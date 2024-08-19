##' get_fit
##' @export
get_fit <- function(sample_name, obj, dipLogR=NA, purity=NA, ploidy=NA, homdel_mb_max=100, neg_mb_max=0, min_ploidy=1.7, max_ploidy=5, min_purity=0.05, max_purity=0.95, cores=1, best_only=T, bin_level=F, purity_stepsize=0.025, ploidy_stepsize=0.025) {
    require(parallel)

    fit_segments <- obj$fit_segments
    sample_dat <- obj$marker_level_annotated[[sample_name]]
    sample_dat <- sample_dat[!is.na(LogR) & !is.na(BAF) & !is.na(LogR_segmented) & !is.na(BAF_segmented)]
    sample_seg <- obj$segment_level[[sample_name]]
    sample_dat <- sample_dat[segment %in% fit_segments,]
    sample_seg <- sample_seg[segment %in% fit_segments,]
    logR_CDF <- obj$ECDF_fits[[sample_name]]$LogR_CDF
    BAF_CDF <- obj$ECDF_fits[[sample_name]]$BAF_CDF

    if(!is.na(dipLogR) & (is.na(purity) | is.na(ploidy))) {        
        message('dipLogR provided')
        auto <- T
        purity <- seq(min_purity,max_purity,by=0.001)
        ploidy <- (2*(purity + 2^(-dipLogR) - 1)) / purity
        myfits <- data.table(pu=purity, pl=round(ploidy,3))
        myfits$dipLogR <- dipLogR
    } else if(!is.na(purity) & !is.na(ploidy)) {
        auto <- F
        message('Purity/ploidy provided')
        myfits <- data.table(pu=purity, pl=ploidy)
        myfits$dipLogR = round(log2(( 2*(1-purity) + purity*2 ) / ( 2*(1-purity) + purity*ploidy )), 4)
    } else if(!is.na(purity) & is.na(ploidy)) {
        message('Purity provided')
        auto <- T
        message('grid-searching ploidy combinations with given purity')
        ploidy <- seq(min_ploidy, max_ploidy, by=ploidy_stepsize)
        myfits <- as.data.table(expand.grid(pu=purity, pl=ploidy))
        myfits$dipLogR = round(log2(( 2*(1-myfits$pu) + myfits$pu*2 ) / ( 2*(1-myfits$pu) + myfits$pu*myfits$pl )), 4)
    } else if(is.na(purity) & !is.na(ploidy)) {
        message('Ploidy provided')
        auto <- T
        message('grid-searching purity combinations with given ploidy')
        purity <- seq(min_purity, max_purity, by=purity_stepsize)
        myfits <- as.data.table(expand.grid(pu=purity, pl=ploidy))
        myfits$dipLogR = round(log2(( 2*(1-myfits$pu) + myfits$pu*2 ) / ( 2*(1-myfits$pu) + myfits$pu*myfits$pl )), 4)
    } else {
        auto <- T
        message('grid-searching purity/ploidy combinations')
        purity <- seq(min_purity, max_purity, by=purity_stepsize)
        ploidy <- seq(min_ploidy, max_ploidy, by=ploidy_stepsize)
        myfits <- as.data.table(expand.grid(pu=purity, pl=ploidy))
        myfits$dipLogR = round(log2(( 2*(1-myfits$pu) + myfits$pu*2 ) / ( 2*(1-myfits$pu) + myfits$pu*myfits$pl )), 4)
    }

    ## for all fits, first filter out impossible fits based on negative copies and too much homozygous deletion
    qc_fit <- function(i, myfits, sample_seg) {
        this.pu <- myfits$pu[i]
        this.pl <- myfits$pl[i]
        this.dipLogR <- myfits$dipLogR[i]
        fit <- get_na_nb(this.pu, this.pl, sample_seg, rfield='LogR_segmented', bfield='BAF_segmented')
        neg_mb <- sum(sample_seg$seg_length[fit$na <= -0.5 | fit$nb <= -0.5],na.rm=T) / 1e6
        homdel_mb <- sum(sample_seg$seg_length[round(fit$na) == 0 & round(fit$nb) == 0],na.rm=T) / 1e6
        na_mb <- sum(sample_seg$seg_length[is.na(fit$na) | is.na(fit$nb)]) / 1e6
        list(pu=this.pu, pl=this.pl, dipLogR=this.dipLogR, neg_mb=neg_mb, homdel_mb=homdel_mb, na_mb=na_mb)
    }
    n_fits <- nrow(myfits)
    profile_list <- mclapply(1:n_fits, qc_fit, myfits, sample_seg, mc.cores=cores) 
    myfits <- rbindlist(profile_list)
    if(auto==T) myfits <- myfits[neg_mb <= neg_mb_max & homdel_mb <= homdel_mb_max & pl <= max_ploidy,]

    get_possible_fits_for_sample <- function(myfits, sample_dat, logR_CDF, BAF_CDF, bin_level, cores) {
        try_fits <- function(i, myfits, sample_dat, logR_CDF, BAF_CDF, bin_level) {
            pu <- myfits$pu[i]
            pl <- myfits$pl[i]
            loglik <- get_loglik(pu, pl, sample_dat, logR_CDF, BAF_CDF, bin_level)
            list(pl=pl, pu=pu, loglik=loglik)
        }
        n_fits <- nrow(myfits)
        l <- mclapply(1:n_fits, try_fits, myfits, sample_dat, logR_CDF, BAF_CDF, bin_level, mc.cores=cores)
        res <- rbindlist(l)
        myfits$loglik <- res$loglik
        myfits
    }

    if(nrow(myfits) > 0) {
        myfits <- get_possible_fits_for_sample(myfits, sample_dat, logR_CDF, BAF_CDF, bin_level, cores)
        myfits <- myfits[order(loglik,decreasing=T),]
        if(best_only==T) myfits <- myfits[1,]
        myfits$sample <- sample_name
        myfits    
    } else {
        NULL
    }
}



##' get_loglik
##' @export
get_loglik <- function(pu, pl, sample_dat, LogR_CDF, BAF_CDF, bin_level=F) {

    get_LogR <- function(pu, pl, na, nb, gc) log2((gc*(1-pu) + pu*(na + nb)) / (gc*(1-pu) + pu*pl))
    get_BAF <- function(pu, pl, na, nb, gc) (1 - pu + pu*nb) / (gc - gc*pu + pu*(na+nb))

    ## get the expected LogR and BAF values if we assume na,nb are actually the nearest integers
    sample_dat <- get_na_nb(p=pu, t=pl, sample_dat, rfield='LogR_segmented', bfield='BAF_segmented')
    #sample_dat <- sample_dat[!is.na(na) & !is.na(nb) & !is.na(BAF) & !is.na(LogR),]
    sample_dat$na <- round(sample_dat$na)
    sample_dat$nb <- round(sample_dat$nb)
    sample_dat[na < 0, na:=0]
    sample_dat[nb < 0, nb:=0]
    sample_dat$LogR_expected <- get_LogR(pu, pl, sample_dat$na, sample_dat$nb, sample_dat$germline_copies)
    sample_dat$BAF_expected <- get_BAF(pu, pl, sample_dat$na, sample_dat$nb, sample_dat$germline_copies)

    if(bin_level==F) {
        ## if doing seg-level likelihood:
        sample_dat <- sample_dat[!duplicated(segment),]
    }
    
    sample_dat <- sample_dat[!is.na(LogR_expected) & !is.na(BAF_expected) & Chromosome %in% c(1:22,'X') & germline_copies > 0,] ## do not use Y or MT for loglik 
    diff_from_LogR_expected <- abs(sample_dat$LogR_segmented - sample_dat$LogR_expected)
    diff_from_LogR_expected <- diff_from_LogR_expected[!is.na(diff_from_LogR_expected)]
    diff_from_BAF_expected <- abs(sample_dat$BAF_segmented - sample_dat$BAF_expected)
    diff_from_BAF_expected <- diff_from_BAF_expected[!is.na(diff_from_BAF_expected)]

    ## calculate the summed loglikelihood of the bins' observed LogRs and BAFs, given the nearest integer copy numbers
    loglik_LogR <- log(c(1-LogR_CDF(diff_from_LogR_expected)))
    min_valid_LogR_loglik <- min(loglik_LogR[!is.infinite(loglik_LogR)])
    max_valid_LogR_loglik <- max(loglik_LogR[!is.infinite(loglik_LogR)])
    loglik_LogR[loglik_LogR==-Inf] <- min_valid_LogR_loglik
    loglik_LogR[loglik_LogR==Inf] <- max_valid_LogR_loglik

    loglik_BAF <- log(c(1-BAF_CDF(diff_from_BAF_expected)))
    min_valid_BAF_loglik <- min(loglik_BAF[!is.infinite(loglik_BAF)])
    max_valid_BAF_loglik <- max(loglik_BAF[!is.infinite(loglik_BAF)])
    loglik_BAF[loglik_BAF==-Inf] <- min_valid_BAF_loglik
    loglik_BAF[loglik_BAF==Inf] <- max_valid_BAF_loglik

    loglik <- sum(loglik_LogR) + sum(loglik_BAF)
    loglik
}


##' plot_fit
##' @export
plot_fit <- function(sample_name, fit, obj, build, int_copies=F, highlight_seg=c(), LogR_min=NA, LogR_max=NA, point_size=0.5) {
    require(cowplot)
    require(ggplot2)
    sample_dat <- obj$marker_level_annotated[[sample_name]]
    sample_seg <- obj$segment_level[[sample_name]]
    sample_seg[, highlight:=F]
    sample_seg[segment %in% highlight_seg, highlight:=T]
    sname <- paste0(sample_name,' (purity=',fit$pu,', ploidy=',fit$pl,', loglik=',round(fit$loglik,3),', dipLogR=',fit$dipLogR,')')
    sex <- obj$ascat.loadData.params$sex
    highlight <- sample_seg[highlight==T,]
    valid_chrs <- unique(sample_seg[germline_copies>0,(Chromosome)])
    build <- obj$ascat.loadData.params$genomeVersion
    chr <- genome_data(build)$chr

    plot_sample_seg <- sample_seg[Chromosome %in% valid_chrs]
    plot_sample_dat <- sample_dat[Chromosome %in% valid_chrs]
    plot_chr <- chr[chr %in% valid_chrs]
    plot_chr$chr <- factor(plot_chr$chr, levels=valid_chrs)
    plot_sample_dat$Chromosome <- factor(plot_sample_dat$Chromosome, levels=valid_chrs)
    plot_sample_seg$Chromosome <- factor(plot_sample_seg$Chromosome, levels=valid_chrs)
    centerline <- median(plot_sample_dat$LogR_segmented,na.rm=T) #median(plot_sample_seg$LogR,na.rm=T) centerline based on segmented median LogR

    LogR_qs <- quantile(plot_sample_dat$LogR, c(0.005,0.995),na.rm=T)
    if(!is.na(LogR_min)) LogR_qs[1] <- LogR_min
    if(!is.na(LogR_max)) LogR_qs[2] <- LogR_max
    plot_sample_dat$arrow <- ''
    plot_sample_dat[LogR > LogR_qs[2], arrow:='up']
    plot_sample_dat[LogR > LogR_qs[2], LogR:=LogR_qs[2]]
    plot_sample_dat[LogR < LogR_qs[1], arrow:='down']
    plot_sample_dat[LogR < LogR_qs[1], LogR:=LogR_qs[1]]
    padding <- (LogR_qs[2] - LogR_qs[1])/50

    p1 <- ggplot(plot_sample_dat, aes(x=global_pos_mb, y=LogR)) + 
        scale_y_continuous(limits=c(LogR_qs[1]-padding, LogR_qs[2]+padding),expand=c(0,0)) + 
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr,expand=c(0,0)) 
    if(nrow(highlight) > 0) p1 <- p1 + annotate("rect", xmin=highlight$global_seg_start_mb, xmax=highlight$global_seg_end_mb, ymin=LogR_qs[1]-padding, ymax=LogR_qs[2]+padding, fill='orange', alpha=0.1)
    p1 <- p1 + geom_point(data=plot_sample_dat[arrow==''], color='#bfbfbf',size=point_size,pch=16) +
        geom_hline(yintercept=centerline, color='green', linewidth=0.25) +
        geom_hline(yintercept=fit$dipLogR, color='purple', linewidth=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end), linewidth=0.25) +
        geom_point(data=plot_sample_dat[arrow=='down'], color='black',fill='#bfbfbf', size=1,pch=25, stroke=0.25) +
        geom_point(data=plot_sample_dat[arrow=='up'], color='black',fill='#bfbfbf', size=1,pch=24, stroke=0.25) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start_mb,xend=global_seg_end_mb,y=LogR_segmented,yend=LogR_segmented),color='blue',linewidth=0.75,lineend='round') +
        theme_fit(base_size=12) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x=NULL, y='LogR', subtitle=sname)

    p2 <- ggplot(plot_sample_dat, aes(x=global_pos_mb, y=BAF)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr,expand=c(0,0)) +
        scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1,by=0.25), expand=c(0,0))
    if(nrow(highlight) > 0) p2 <- p2 + annotate("rect", xmin=highlight$global_seg_start_mb, xmax=highlight$global_seg_end_mb, ymin=-0.05, ymax=1.05, fill='orange', alpha=0.1)
    p2 <- p2 + geom_point(color='#bfbfbf',size=point_size,pch=16) +
        geom_vline(xintercept=c(0,plot_chr$global_end), linewidth=0.25) +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start_mb,xend=global_seg_end_mb,y=BAF_segmented,yend=BAF_segmented),color='blue',linewidth=0.75,lineend='round') +
        geom_segment(data=plot_sample_seg,aes(x=global_seg_start_mb,xend=global_seg_end_mb,y=1-BAF_segmented,yend=1-BAF_segmented),color='blue',linewidth=0.75,lineend='round') +
        theme_fit(base_size=12) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x=NULL, y='BAF')

    ascn <- get_na_nb(p=fit$pu, t=fit$pl, x=plot_sample_seg, rfield='LogR_segmented', bfield='BAF_segmented')
    ascn[,tcn:=na+nb]
    ascn[na >= nb,mcn:=nb]
    ascn[na < nb,mcn:=na]
    ascn[is.na(nb) & !is.na(na), tcn:=na]
    ascn[is.na(nb) & !is.na(na), mcn:=NA]
    ascn[tcn >= 10, tcn:=log10(tcn)+9]
    ascn[mcn >= 10, tcn:=log10(mcn)+9]
    ascn[,tcn_ad_from_int:=abs(tcn-round(tcn))] 
    ascn[,mcn_ad_from_int:=abs(mcn-round(mcn))] 
    ascn[tcn < -1, tcn:=-1]
    ascn[mcn < -1, mcn:=-1]

    if(int_copies==T) {
        ascn[,tcn:=round(tcn)]
        ascn[,mcn:=round(mcn)]
    }

    p3 <- ggplot(ascn) + 
        scale_y_continuous(limits=c(-1,11.5), breaks=c(seq(0,10,by=2),11), labels=c(seq(0,10,by=2),100), expand=c(0,0)) +
        scale_x_continuous(breaks=plot_chr$global_midpoint, labels=plot_chr$chr, expand=c(0,0)) 
    if(nrow(highlight) > 0) p3 <- p3 + annotate("rect", xmin=highlight$global_seg_start_mb, xmax=highlight$global_seg_end_mb, ymin=-Inf, ymax=Inf, fill='orange', alpha=0.1)
    p3 <- p3 + 
        geom_hline(yintercept=seq(0,11,by=1), color='#bfbfbf', linewidth=0.25) +
        geom_vline(xintercept=c(0,plot_chr$global_end), linewidth=0.25) +
        geom_segment(aes(x=global_seg_start_mb,xend=global_seg_end_mb,y=tcn,yend=tcn),linewidth=1,color='black',lineend='round') +
        geom_segment(aes(x=global_seg_start_mb,xend=global_seg_end_mb,y=mcn,yend=mcn),linewidth=1,color='red',lineend='round') +
        theme_fit(base_size=12) +
        theme(legend.position='bottom') +
        labs(x='Genomic posititon', y='Copy number') 
    p <- plot_grid(p1, p2, p3, align='v', ncol=1, rel_heights=c(1.05,1,1.2), axis='lr')
    p
}


get_na_nb <- function(p, t, x, rfield='LogR', bfield='BAF', gfield='germline_copies') {
    r <- x[[rfield]]
    b <- x[[bfield]]
    g <- x[[gfield]]
    c <- 2^r
    c1 <- 2^(r+1)

    ## NB: here we model the median admixed ploidy as 2*(1-p)+p*t
    ## this works because the LogR values are calculated as log2(x/median(x))
    na = -1*( b*p*c*t - b*p*c1 + b*c1 - g*p + g - p*c*t + p*c1 + p - c1 - 1 ) / p
    nb = (b*p*c*t - b*p*c1 + b*c1 + p - 1) / p        

    ## go over previous values where BAF was NA
    NAs <- is.na(b)
    tmp_na <- (g*(p-1) + c*(p*(t-2)+2)) / p
    na[NAs] <- tmp_na[NAs]
    nb[NAs] <- as.numeric(NA)
    x$na <- na
    x$nb <- nb
    x
}


##' theme_fit
##' @export
theme_fit <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size,
        base_rect_size = base_rect_size) %+replace% theme(line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"),
        text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = ggplot2::margin(), debug = F), 
        axis.text = element_text(colour = "black", size = rel(0.8)), 
        axis.ticks = element_line(colour = "black", linewidth = rel(1)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(), #line(colour = "black", size = rel(1)), 
        legend.key = element_blank()) 
} 


