##' plot_tcn_angular_distance_heatmap
##'
##' Plot a preliminary copy number heatmap based on angular distances from total copy number
##'
##' @export
plot_tcn_angular_distance_heatmap <- function(obj=NULL, segs=NULL, build=NULL, sample_groups=NULL, group_colors=NULL, normal_sample, patient='', squish_tree_factor=1.25) {
    
    if(!is.null(obj)) {
        segs <- rbindlist(obj$segment_level, fill=T)
        build <- obj$main_params$build
        patient <- obj$main_params$patient
        patient <- obj$main_params$patient
    } 
    
    sex <- obj$main_params$sex
    if(sex=='XX') {
        valid_chr <- c(1:22,'X')
    } else {
        valid_chr <- c(1:22,'X','Y')
    }
    
    if(is.null(obj) & (is.null(segs) | is.null(build))) stop('Must provide either [obj] OR ([segs] AND [build])') 
    
    ## angular distance tree from segment LogRs
    LogR_mat <- d2m(data.table::dcast(segment ~ sample, value.var='LogR_segmented', data=segs))
    LogR_mat <- LogR_mat[!is.na(rowSums(LogR_mat)),]
    problematic_samples <- colnames(LogR_mat)[is.na(colSums(LogR_mat))]
    if(length(problematic_samples) > 0) {
        message('Sample with NA LogR values detected!')
        message(paste(problematic_samples, collapse=', '))
        message('Removing it and proceeding.')
    }
    LogR_mat <- LogR_mat[,!colnames(LogR_mat) %in% problematic_samples]
    
    n_seg <- nrow(LogR_mat)
    n_samp <- ncol(LogR_mat)
    for(i in 1:n_samp) LogR_mat[,i] <- LogR_mat[,i] / sqrt(sum(LogR_mat[,i]^2))
    
    long <- as.data.table(reshape2::melt(LogR_mat))
    names(long) <- c('segment','sample','value')
    segs <- merge(long, obj$segments, by='segment', all.x=T)
    segs <- segs[Chromosome %in% valid_chr,]
    segs$Chromosome <- factor(segs$Chromosome, levels=valid_chr)
    
    ## get angular distance matrix/tree
    numerator <- t(LogR_mat) %*% LogR_mat
    magnitudes <- sqrt(colSums(LogR_mat^2))
    denominator <- magnitudes %*% t(magnitudes)
    cos_theta_LogR <- numerator / denominator
    theta_LogR <- acos(cos_theta_LogR)
    for(i in 1:ncol(theta_LogR)) theta_LogR[i,i] <- 0
    theta_LogR <- cbind(theta_LogR, diploid=pi/3)
    theta_LogR <- rbind(theta_LogR, diploid=pi/3)
    for(i in 1:(ncol(LogR_mat)+1)) theta_LogR[i,i] <- 0
    tree_LogR <- nj(theta_LogR)
    tree_LogR <- phytools::reroot(tree_LogR, which(tree_LogR$tip.label=='diploid'))
    tree_LogR$tip.label[tree_LogR$tip.label=='diploid'] <- normal_sample
    
    setnames(segs,c('Chromosome'),c('chr'))
    segs$sample <- as.character(segs$sample)
    message('Expanding segments to include NA regions in each chromosome ...')
    
    expand_segments_to_complete_chromosome_for_sample <- function(this.sample, segs, gr_chr) {
        message(this.sample)
        mat_sample <- segs[sample==this.sample,]
        mat_sample$chr <- factor(mat_sample$chr, levels(seqnames(gr_chr)))
        mat_sample <- mat_sample[order(chr)]
        gr_mat_sample <- makeGRangesFromDataFrame(mat_sample,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='seg_start', end.field='seg_end')
        NA_regions <- BiocGenerics::setdiff(gr_chr, gr_mat_sample)
        complete_regions <- as.data.table(sort(c(gr_mat_sample, NA_regions)))
        complete_regions$sample <- this.sample
        complete_regions[,segment:=paste0(seqnames,':',start,'-',end)]
        setnames(complete_regions,'seqnames','chr')
        complete_regions <- complete_regions[order(chr, start, end),]
        ## collapse regions with no difference in copy number
        complete_regions$sample <- this.sample
        complete_regions
    }
    
    gd <- obj$gd$chr
    gd <- gd[chr %in% valid_chr]
    gd$chr <- factor(gd$chr, levels=valid_chr)
    gd[,chr_start:=chr_start * 1e6]
    gd[,chr_end:=chr_end * 1e6]
    
    gr_chr <- makeGRangesFromDataFrame(gd,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='chr_start', end.field='chr_end')
    sample_list <- lapply(unique(segs$sample), expand_segments_to_complete_chromosome_for_sample, segs, gr_chr)
    segs <- rbindlist(sample_list)
    
    segs <- merge(segs, gd[,c('chr','global_start'),with=F], by='chr', all.x=T)
    segs[,global_seg_start_mb:=global_start + start/1e6]
    segs[,global_seg_end_mb:=global_start + end/1e6]
    segs <- segs[chr %in% valid_chr]
    segs$chr <- factor(segs$chr, levels=valid_chr)
    segs[,segment:=paste0(chr,':',start,'-',end)]
    
    ## add values for the diploid normal
    toadd <- segs[!duplicated(segment),]
    toadd[,value:=0]
    toadd[,sample:=normal_sample]
    segs <- rbind(segs, toadd)
    
    p1 <- ggtree(tree_LogR, linewidth=0.5)
    p1 <- p1 + guides(fill='none') + theme(legend.position='none') + labs(title=patient, subtitle='CN angular distance tree')
    pd <- as.data.table(p1$data)
    pd <- pd[isTip==T,]
    pd <- pd[order(y,decreasing=F),]
    pd$label <- factor(pd$label, levels=pd$label)
    segs$sample <- factor(segs$sample, levels=pd$label)
    segs$sample.i <- as.integer(segs$sample)
    
    if(!is.null(sample_groups)) {
        message('Adding sample type sample_groups')
        p1 <- p1 %<+% sample_groups[,c('lesion_id','group'),with=F]
        p1 <- p1 + geom_tiplab(aes(color=group))
        p1 <- p1 + scale_color_manual(values=group_colors,name='Organ category')
    } else {
        p1 <- p1 + geom_tiplab()
    }
    
    p1 <- p1 + xlim(0, squish_tree_factor*max(p1$data$x))
    gd2 <- copy(gd)
    chr19_22 <- gd[chr==19,]
    chr19_22$length <- sum(gd[chr %in% c(19,20,21,22),(length)])
    chr19_22[,global_midpoint:=global_start+0.5*length]
    chr19_22[,global_end:=global_start+length]
    chr19_22$chr <- '19-22'
    gd2 <- gd2[!chr %in% c(19,20,21,22)]
    gd2 <- rbind(gd2, chr19_22) 
    gd2 <- gd2[order(global_midpoint),]
    
    right_labs <- data.table(sample=levels(segs$sample))
    right_labs$sample <- factor(right_labs$sample, levels=right_labs$sample)
    right_labs$sample.i <- as.integer(right_labs$sample)
    right_labs[,y:=sample.i]
    
    p2 <- ggplot(segs) + 
        scale_x_continuous(expand=c(0,-0.1), breaks=gd2$global_midpoint, labels=gd2$chr) + 
        scale_y_continuous(expand=c(0,0), breaks=right_labs$y, labels=right_labs$sample, position='right') + 
        geom_rect(aes(ymin=sample.i-0.5, ymax=sample.i+0.5, xmin=global_seg_start_mb, xmax=global_seg_end_mb, fill=value)) + 
        scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, na.value='#bfbfbf',name='Normalized TCN') +
        geom_vline(xintercept=c(0,tail(gd$global_end,1)), linewidth=0.5, linetype='solid') +
        geom_hline(yintercept=seq(min(segs$sample.i)-0.5,max(segs$sample.i)+0.5), linewidth=0.25, linetype='solid') + 
        geom_vline(xintercept=tail(gd$global_start,-1), linewidth=0.25, linetype='dotted') + 
        guides(color='none') +
        theme_ang(base_size=10) +
        labs(x='Genomic Position', y=NULL, title='', subtitle='Normalized total copy number') +
        theme(axis.line=element_blank(), legend.position='bottom', axis.ticks=element_blank(), axis.text.y=element_blank())
    
    p <- plot_grid(p1, p2, nrow=1, rel_widths=c(1.5,3), align='h', axis='tb')
    p
}
