
##' preprocess_bin_data
##'
##' Create preprocessed data from snp-pileup, QDNAseq, and GLIMPSE output.
##'
##' @export
preprocess_bin_data <- function(qdnaseq_data, pileup_data, phased_bcf, sample_map, normal_sample, sex, build, max_phaseable_distance, min_bin_reads_for_baf, blacklisted_regions_file, LogR_range_allowed, LogR_winsor_percentiles, LogR_smooth_bins, normal_correction, imbalance_alpha=0.05) { 

    all_chrs <- c(1:22,'X','Y','MT')
    if(sex=='XX') {
        diploid_chrs <- c(1:22,'X')
    } else {
        diploid_chrs <- c(1:22)       
    }
    
    ## load bin counts from QDNAseq. We need to annotate these bins with allele-specific counts (for bins containing het-SNPs)
    cnt <- readRDS(qdnaseq_data)
    setnames(cnt,c('start','end'),c('bin_start','bin_end'))
    cnt[,bin_length:=bin_end-bin_start+1]

    if(!is.na(blacklisted_regions_file)) {
        message('Excluding blacklisted regions and correcting bin counts ...')
        ## scale the counts based on the length of removed targeted regions
        targets <- fread(blacklisted_regions_file)
        targets <- targets[,c(1:4)]
        names(targets) <- c('chr','target_start','target_end','target_id')
        setkey(cnt,'chromosome','bin_start','bin_end') 
        setkey(targets,'chr','target_start','target_end') 
        cnt <- foverlaps(cnt, targets, type='any')
        cnt[,target_length:=target_end-target_start]
        cnt[is.na(target_length), target_length:=0]

        ## some bins have multiple targets, collapse them and 
        ## subtract their combined length from the bin length
        collapse_cnt_targets <- function(cnt) {
            total_target_length <- sum(cnt$target_length)
            total_bin_length <- unique(cnt$bin_length)
            targets <- cnt$target_id[!is.na(cnt$target_id)]
            targets <- paste(targets, collapse=', ')
            list(bin_length=total_bin_length - total_target_length, targets=targets)
        }
        cnt <- cnt[,collapse_cnt_targets(.SD), by=c('chromosome','feature','bin_start','bin_end','count','file')]    
    }

    ## normalize the counts for bin-lenth (inluding after removing any blacklisted regions)
    mu_length <- mean(cnt$bin_length)
    cnt$count_per_bp <- cnt$count / cnt$bin_length
    cnt$count <- cnt$count_per_bp * mu_length

    ## update sample names
    map <- fread(sample_map)
    cnt <- merge(cnt, map[,c('qdnaseq_name','proper_name'),with=F], by.x='file', by.y='qdnaseq_name', all.x=T)
    cnt <- cnt[!is.na(proper_name)]
    cnt[,file:=NULL]
    setnames(cnt,c('feature','chromosome','proper_name'),c('bin','chr','sample'))

    ## format and order the data
    cnt$chr <- factor(cnt$chr, levels=all_chrs)
    cnt <- cnt[!is.na(chr),]
    cnt <- cnt[order(chr,bin_start,bin_end),]

    ## get the names of the samples
    samples <- map$proper_name
    tumor_samples <- samples[!samples %in% normal_sample]
    
    ## load phased SNPs
    message('Loading phased SNP bcf-file: ',phased_bcf,' ...')
    bcf <- scanBcf(phased_bcf)
    phased_snps <- data.table(CHROM=bcf$CHROM, POS=bcf$POS, ID=bcf$ID, REF=bcf$REF, ALT=bcf$ALT, PHASE=bcf$GENO$GT[,1])

    ## subset for chrs that are diploid in normal given sex
    phased_snps$CHROM <- gsub('chr','',phased_snps$CHROM)
    phased_snps <- phased_snps[CHROM %in% diploid_chrs,]
    phased_snps[,POS2:=POS+1]
    phased_snps$CHROM <- factor(phased_snps$CHROM, levels=diploid_chrs) ## use all chrs for levels
    setkey(phased_snps,'CHROM','POS','POS2')

    ## annotate the phased snps with binned counts
    bins <- cnt[!duplicated(bin),c('bin','chr','bin_start','bin_end'),with=F]
    setkey(bins,'chr','bin_start','bin_end')

    phased_snps <- foverlaps(phased_snps, bins, type='within')
    phased_snps <- phased_snps[!is.na(bin),]
    phased_snps$bin <- factor(phased_snps$bin, levels=unique(phased_snps$bin))
    phased_snps$bin.i <- as.integer(phased_snps$bin)

    ## define phasing blocks
    message('Defining phasing blocks with maximum phaseable distance: ',max_phaseable_distance,' ...')
    phased_snps <- block_phased_snps(phased_snps, max_phaseable_distance=max_phaseable_distance)

    ## get the number of SNPs per phasing block in each bin
    collapse_blocks <- function(phased_snps) {
        block_start = min(phased_snps$POS)
        block_end = max(phased_snps$POS)
        block_snps = length(unique(phased_snps$POS))
        list(block_start=block_start, block_end=block_end, block_snps=block_snps)
    }
    blocks <- phased_snps[, collapse_blocks(.SD), by=c('CHROM','bin','block')]
    blocks[,block_midpoint:=round((block_start+block_end)/2)]

    ## load allele-specific read counts for the phased het-SNPs for each sample
    message('Loading snp-pileup data: ',pileup_data,' ...')
    d <- fread(pileup_data, sep=',')
    d$Chromosome <- gsub('chr','',d$Chromosome)
    d <- d[Chromosome %in% diploid_chrs]
    d$Chromosome <- factor(d$Chromosome, levels=diploid_chrs)
    d <- merge(d, phased_snps[,c('CHROM','POS','PHASE','bin','block'),with=F], by.x=c('Chromosome','Position'), by.y=c('CHROM','POS'), all.x=T)
    d <- d[!is.na(block),]
    front_fields <- names(d)[grepl('File',names(d))==F]
    back_fields <- names(d)[!names(d) %in% front_fields]
    d <- d[,c(front_fields, back_fields),with=F]

    ## replace the 'File' names with sample-names
    tmp <- data.table(pos=1:ncol(d), orig_name=names(d))
    tmp[grepl('File', orig_name),pileup_pos:=as.integer(gsub('File','',strtrim(orig_name,nchar(orig_name)-1)))]
    tmp[grepl('File', orig_name),allele:=substr(orig_name,nchar(orig_name), nchar(orig_name))]
    tmp <- merge(tmp, map[,c('pileup_pos','proper_name'),with=F], by='pileup_pos', all.x=T)
    tmp <- tmp[order(pos),]
    tmp[!is.na(pileup_pos), new_name:=paste0(proper_name,'.',allele)] 
    tmp[is.na(pileup_pos), new_name:=orig_name]
    names(d) <- tmp$new_name 

    ## extract and combine read counts per block for each allele
    count10_A <- d[PHASE=='1|0',c('bin','block',paste0(samples,'.A')),with=F]; names(count10_A) <- c('bin','block',samples)
    count01_R <- d[PHASE=='0|1',c('bin','block',paste0(samples,'.R')),with=F]; names(count01_R) <- c('bin','block',samples)
    count1 <- rbind(count10_A, count01_R)
    count01_A <- d[PHASE=='0|1',c('bin','block',paste0(samples,'.A')),with=F]; names(count01_A) <- c('bin','block',samples)
    count10_R <- d[PHASE=='1|0',c('bin','block',paste0(samples,'.R')),with=F]; names(count10_R) <- c('bin','block',samples)
    count2 <- rbind(count01_A, count10_R)
    collapse_block <- function(count) {
        out <- colSums(count)     
        as.list(out)
    }
    block_counts1 <- count1[,collapse_block(.SD),by=c('bin','block')]
    block_counts2 <- count2[,collapse_block(.SD),by=c('bin','block')]
    block_counts <- merge(block_counts1, block_counts2, by=c('bin','block'), all=T)
    block_counts[is.na(block_counts)] <- 0

    x_samples <- paste0(samples,'.x')
    y_samples <- paste0(samples,'.y')
    x_counts <- block_counts[,c('bin','block',x_samples),with=F]
    y_counts <- block_counts[,c('bin','block',y_samples),with=F]
    x_counts <- data.table::melt(x_counts, id.vars=c('bin','block'))
    x_counts[,variable:=gsub('[.]x$','',as.character(variable))]
    y_counts <- data.table::melt(y_counts, id.vars=c('bin','block'))
    y_counts[,variable:=gsub('[.]y$','',as.character(variable))]
    block_counts_long <- merge(x_counts, y_counts, by=c('bin','block','variable'))
    block_counts_long[,total:=value.x+value.y]
    
    ## test each block within a bin for significant allelic imbalance.
    ## If any do, we assume there is the same allelic imbalance across all the blocks.
    ## In this case, we will swap .x and .y labels for the other blocks so that their combined BAFs have the same direction.
    test_blocks <- function(qc,tumor_samples) { 
        x <- as.numeric(qc[,paste0(tumor_samples,'.x'),with=F]) 
        y <- as.numeric(qc[,paste0(tumor_samples,'.y'),with=F]) 
        z <- y - x
        n <- length(z)
        suppressWarnings(w <- wilcox.test(z, mu = 0, exact = FALSE, correct = FALSE))
        p = as.numeric(w$p.value)
        W = w$statistic
        mean_W = n * (n + 1) / 4
        sd_W = sqrt(n * (n + 1) * (2 * n + 1) / 24)
        z_score = (W - mean_W) / sd_W
        effsize = z_score / sqrt(n)
        dp <- sum(x) + sum(y)
        list(p=p, effsize=effsize, dp=dp)
    }
    block_summary <- block_counts[,test_blocks(.SD, tumor_samples),by=c('bin','block')]

    get_block_fdr_per_bin <- function(block_summary, alpha) {
        ## annotate all blocks in a bin with the q-val and LFC for the MOST SIGNIFICANT block
        ## if there are multiple significant blocks with LFCs in different directions, note this as problematic
        ## if there were no significant blocks in the bin, we will discard the BAF for this bin
        block_summary$q <- p.adjust(block_summary$p, method='BH')
        n_sig_neg <- sum(block_summary$q < alpha & block_summary$effsize < 0,na.rm=T)
        n_sig_pos <- sum(block_summary$q < alpha & block_summary$effsize > 0,na.rm=T)
        out <- block_summary[order(q, decreasing=F),]
        if(n_sig_neg + n_sig_pos == 0) {
            out$result <- 'no significant imbalance'
        } else if(n_sig_neg > 0 & n_sig_pos > 0) {
            out$result <- 'discordant signal'
        } else {
            out$result <- 'imbalanced'
        }
        out[,c('p','q','effsize','result','dp')]
    }
    block_summary <- block_summary[,get_block_fdr_per_bin(.SD, alpha=imbalance_alpha),by=block]

    ## in every bin with a significant block, annotate its q-val and effsize (NB: effsize > 0 means .y > .x)
    block_counts <- merge(block_counts, block_summary, by='block', all.x=T)

    # discard the bins with discordant signals
    block_counts <- block_counts[result!='discordant signal',]
    
    # in bins with multiple blocks, choose the most-significant as the representative block
    collapsed_bin_counts = block_counts[order(bin, q, decreasing=F),]
    collapsed_bin_counts <- collapsed_bin_counts[!duplicated(bin),]
    collapsed_bin_counts <- collapsed_bin_counts[order(bin),]
    setnames(collapsed_bin_counts,c('p','q','effsize','result','dp'),c('bin_p','bin_q','bin_effsize','bin_status','bin_dp')) 

    ## get BAFs for each bin/sample based on x,y counts
    tumors.x <- paste0(tumor_samples,'.x')
    tumors.y <- paste0(tumor_samples,'.y')
    normal.x <- paste0(normal_sample,'.x')
    normal.y <- paste0(normal_sample,'.y')
    samples.x <- paste0(samples,'.x')
    samples.y <- paste0(samples,'.y')
    
    res.x <- melt(collapsed_bin_counts[,c('bin','bin_p','bin_q','bin_effsize','bin_status','bin_dp',samples.x),with=F], id.vars=c('bin','bin_p','bin_q','bin_effsize','bin_status','bin_dp'))
    res.x[,variable:=gsub('[.]x','',as.character(variable))] 
    setnames(res.x,'value','count.x')
    res.y <- melt(collapsed_bin_counts[,c('bin',samples.y),with=F], id.vars=c('bin'))
    res.y[,variable:=gsub('[.]y','',as.character(variable))] 
    setnames(res.y,'value','count.y')
    res <- merge(res.x, res.y, by=c('bin','variable'), all=T) 
    res[,dp:=count.x+count.y]
    setnames(res,'variable','sample')
    res[,baf:=count.y/dp] 

    ## combine collapsed bin bafs with the count data
    d <- merge(cnt, res, by=c('sample','bin'), all.x=T)
    d <- d[order(chr,bin_start,bin_end),]
    d$bin <- factor(d$bin, levels=unique(d$bin)) 
    d[dp < min_bin_reads_for_baf, baf:=NA]
    dN <- d[sample==normal_sample,c('bin','count')]
    setnames(dN,'count','n_count')
    d <- merge(d, dN, by=c('bin'), all.x=T)
    d[,position:=round((bin_start + bin_end - 1)/2)]
    d <- d[order(chr,position,sample),]

    message('Adding chr-arms ...')
    arms <- genome_data(build)$arms
    arms[,arm_start:=arm_start*1e6]  
    arms[,arm_end:=arm_end*1e6]  
    setkey(arms,'chr','arm_start','arm_end')
    d[,start:=(position/1e6) - 0.5]
    d[,end:=(position/1e6) + 0.5]
    setkey(d,'chr','start','end')
    d <- foverlaps(d, arms, type='within')
    setnames(d,c('chr','position'),c('Chromosome','Position'))
    d[Chromosome=='MT', arm:='']
    d <- d[!is.na(arm),]
    d[,c('start','end','arm_start','arm_end'):=NULL]
    d[,charm:=paste0(Chromosome,arm)]

    if(normal_correction==0) {
        ## calculate LogR, first by the normal depth, then by the sample median T/N, but only among the autosomes
        .get_LogR <- function(d) {
            mid <- median(d$count[d$Chromosome %in% c(1:22)],na.rm=T)
            d$LogR <- log2(d$count / mid)
            d
        }
        d <- d[,.get_LogR(.SD),by=sample]

    } else {
        ## calculate LogR, first by the normal depth, then by the sample median T/N, but only among the autosomes
        .get_LogR <- function(d) {
            Ratio <- d$count / d$n_count
            mid <- median(Ratio[d$Chromosome %in% c(1:22)],na.rm=T)
            d$LogR <- log2(Ratio / mid)
            d
        }
        d <- d[,.get_LogR(.SD),by=sample]

        # from ASCAT https://github.com/VanLoo-lab/ascat/blob/master/ASCAT/R/ascat.prepareHTS.R 
        # For males, chrX needs to be adjusted as logR baseline will be 0 because of T/N ratio
        if (sex=="XY") {

            message('Adjusting LogR in chrX,Y PARs for male.')
            # PAR1 and PAR2 information should be a mix of chrX and chrY so we should expect 1+1 (1 copy from X and 1 copy from Y).
            # nonPAR should be X-specific and baseline is 1+0 so logR needs to be decreased according to gamma parameter (ascat.runAscat)
            if (build=="hg19") {
                XnonPAR=c(2699521, 154931043)
                YnonPAR=c(2699521, 59034049)
            } else if (build=="hg38") {
                XnonPAR=c(2781480, 155701382)
                YnonPAR=c(2781480, 56887902)
            }
            XnonPAR=which(d$Chromosome %in% c("X", "chrX") & d$Position>=XnonPAR[1] & d$Position<=XnonPAR[2])
            YnonPAR=which(d$Chromosome %in% c("Y", "chrY") & d$Position>=YnonPAR[1] & d$Position<=YnonPAR[2])
            d$LogR[XnonPAR]=d$LogR[XnonPAR]-1
            d$LogR[YnonPAR]=d$LogR[YnonPAR]-1
        }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # These options can help accommodate both 
    # off-target read data and true lpWGS data
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(!is.na(LogR_range_allowed[1])) {
        message(paste0('Removing LogR values < ',LogR_range_allowed[1],' ...'))
        d[Chromosome!='MT' & LogR < LogR_range_allowed[1], LogR:=NA]
    }
    
    if(!is.na(LogR_range_allowed[2])) {
        message(paste0('Removing LogR > ',LogR_range_allowed[2],' ...'))
        d[Chromosome!='MT' & LogR > LogR_range_allowed[2], LogR:=NA]
    }

    if(!is.na(LogR_winsor_percentiles[1])) {
        q <- quantile(d$LogR[d$Chromosome!='MT'],LogR_winsor_percentiles[1],na.rm=T) 
        message(paste0('Winsorising LogR below ',LogR_winsor_percentiles[1],' percentile (LogR < ',q,') ...'))
        d[Chromosome!='MT' & LogR < q, LogR:=q]
    }
    
    if(!is.na(LogR_winsor_percentiles[2])) {
        q <- quantile(d$LogR[d$Chromosome!='MT'],LogR_winsor_percentiles[2],na.rm=T) 
        message(paste0('Winsorising LogR above ',LogR_winsor_percentiles[2],' percentile (LogR > ',q,') ...'))
        d[Chromosome!='MT' & LogR > q, LogR:=q]
    }

    if(!is.na(LogR_smooth_bins)) {
        message('Smoothing LogR with a symmetric rolling median (+/- ',LogR_smooth_bins,' bins) ...')
        get_smooth_LogR <- function(d,width) {
            n <- nrow(d)
            d <- d[order(bin_start,bin_end)]
            d$LogRsmooth <- as.numeric(NA)
            for(i in 1:n) {
                indices <- (i-width):(i+width)
                indices <- indices[indices >= 1 & indices <= n]
                rolling_median <- median(d$LogR[indices],na.rm=T)
                d$LogRsmooth[i] <- rolling_median
            }
            d
        }
        dm <- d[Chromosome=='MT',]
        d <- d[Chromosome!='MT',get_smooth_LogR(.SD,width=LogR_smooth_bins),by=c('sample','Chromosome','arm')]
        d[,LogR:=LogRsmooth]
        d[,LogRsmooth:=NULL]
        d <- rbind(d, dm)
        message('Done with smoothing.')
    }
    
    ## add count1/2, which are really the normalized counts from QDNAseq scaled by BAF and (1-BAF) for the respective bin
    d[is.nan(baf),baf:=NA]
    d[!is.na(baf),count1:=count * baf]
    d[!is.na(baf),count2:=count * (1-baf)]
    d[is.na(baf),count1:=count]
    d[is.na(baf),count2:=0]
    d[is.na(baf),Ref:='.']
    d[is.na(baf),Alt:='.']
    d[!is.na(baf),Ref:='N']
    d[!is.na(baf),Alt:='N']
    setnames(d,'baf','BAF')
    d[,dp:=count1+count2]
    d <- d[,c('Chromosome','arm','Position','Ref','Alt','sample','count1','count2','LogR','BAF','bin','charm'),with=F]
    d[BAF < 0 | BAF > 1, BAF:=NA]
    d[,bin:=paste0(Chromosome,':',Position)]

    ## make sure all bin positions are valid
    chr <- genome_data(build)$chr
    d <- merge(d, chr[,c('chr','global_start','global_end'),with=F], by.x='Chromosome', by.y='chr', all.x=T)
    d[,global_pos:=global_start + (Position/1e6)]
    d <- d[global_pos >= global_start & global_pos <= global_end,]
    d$Chromosome <- factor(d$Chromosome, levels=all_chrs)
    d <- d[order(Chromosome,Position),]
    d$bin <- as.integer(factor(d$bin, levels=unique(d$bin)))
    message('Done!')
    d
}
