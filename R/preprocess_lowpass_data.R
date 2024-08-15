##' preprocess_lowpass_data
##' Create preprocessed data from output snp-pileup, QDNAseq, and GLIMPSE
##' @export
preprocess_lowpass_data <- function(qdnaseq_data, pileup_data, phased_bcf, sample_map, sex, build='hg19', max_phaseable_distance=20000, min_block_reads=10, seed=NA, normal_sample, min_tumors_for_imputing=NA, blacklisted_regions_file=NA) {

    if(!is.na(seed)) set.seed(seed)

    all_chrs <- c(1:22,'X','Y','MT')
    if(sex=='XX') {
        diploid_chrs <- c(1:22,'X')
    } else {
        diploid_chrs <- c(1:22)       
    }
    
    ## load bin counts from QDNAseq. We need to annotate these bins with allele-specific counts 
    ## (for bins containing het-SNPs)
    ## this probably means the het-SNPs bins need to perfectly align with the QDNAseq bins
    cnt <- readRDS(qdnaseq_data)
    setnames(cnt,c('start','end'),c('bin_start','bin_end'))
    cnt[,bin_length:=bin_end-bin_start+1]

    if(!is.na(blacklisted_regions_file)) {
        message('Excluding blacklisted regions and correcting bin counts ...')
        ## scale the counts based on the length of removed targeted regions
        targets <- fread(blacklisted_regions_file)
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

    ## test each block within a bin for significant allelic imbalance.
    ## If any do, we assume there is the same allelic imbalance across all the blocks.
    ## In this case, we will swap .x and .y labels for the other blocks so that their combined BAFs have the same direction.
    test_blocks <- function(qc,tumor_samples) { 
        x <- as.numeric(qc[,paste0(tumor_samples,'.x'),with=F]) 
        y <- as.numeric(qc[,paste0(tumor_samples,'.y'),with=F]) 
        sum_x <- sum(x)
        sum_y <- sum(y)
        suppressWarnings(p <- wilcox.test(x, y, paired=T, alternative='two.sided')$p.value)
        lfc <- log2(mean(y) / mean(x))
        list(p=p, lfc=lfc, sum_x=sum_x, sum_y=sum_y)
    }
    block_summary <- block_counts[,test_blocks(.SD, tumor_samples),by=c('bin','block')]

    correct_blocks <- function(block_summary, alpha=0.1) {
        ## annotate all blocks in a bin with the q-val and LFC for the MOST SIGNIFICANT block
        ## if there are multiple significant blocks with LFCs in different directions, note this as problematic
        block_summary$q <- p.adjust(block_summary$p, method='BH')
        block_summary$problematic <- any(block_summary[q < alpha,(lfc)] < 0) & any(block_summary[q < alpha,(lfc) > 0])
        block_summary <- block_summary[order(q, decreasing=F),]
        block_summary[1,c('p','q','lfc','problematic')]
    }
    block_summary <- block_summary[,correct_blocks(.SD),by=bin]
    block_summary[problematic==T,lfc:=NaN]
    block_summary[problematic==T,q:=NaN]

    ## in every bin with a significant block, annotate its q-val and LFC (NB: LFC > 0 means .y > .x)
    setnames(block_summary,c('p','q','lfc'),c('bin_p','bin_q','bin_lfc'))
    block_counts <- merge(block_counts, block_summary[,c('bin','bin_p','bin_q','bin_lfc'),with=F], by='bin', all.x=T)
    x_samples <- paste0(samples,'.x')
    y_samples <- paste0(samples,'.y')
    x_counts <- block_counts[,(x_samples),with=F]
    y_counts <- block_counts[,(y_samples),with=F]
    block_counts$block_lfc <- log2(rowSums(y_counts) / rowSums(x_counts))
    block_counts[is.nan(bin_lfc), bin_lfc:=NA]
    block_counts[is.nan(block_lfc), block_lfc:=NA]

    ## reorient blocks so that their x,y orientations match that of the most-significant block in each bin 
    reorient_blocks <- function(current_block_counts, x_samples, y_samples) {
        different_directions <- current_block_counts$block_lfc * current_block_counts$bin_lfc < 0
        if(different_directions==T & !is.na(different_directions)) { 
            current_x_counts <- current_block_counts[,(x_samples),with=F]
            current_y_counts <- current_block_counts[,(y_samples),with=F]
            current_block_counts[,(x_samples)] <- current_y_counts
            current_block_counts[,(y_samples)] <- current_x_counts
        }
        current_block_counts       
    } 
    realigned_block_counts <- block_counts[,reorient_blocks(.SD, x_samples, y_samples),by=block]
    new_x_counts <- realigned_block_counts[,(x_samples),with=F]
    new_y_counts <- realigned_block_counts[,(y_samples),with=F]
    realigned_block_counts$realigned_block_lfc <- log2(rowSums(new_y_counts) / rowSums(new_x_counts))

    ## collapse to bin-level data. For now, just combine the .x and .y read counts across all blocks regardless of bin-level significance
    collapse_bins <- function(realigned_block_counts, x_samples, y_samples) {
        out <- realigned_block_counts[1,c('bin_p','bin_q','bin_lfc'),with=F]
        current_x_counts <- realigned_block_counts[,(x_samples),with=F]
        current_y_counts <- realigned_block_counts[,(y_samples),with=F]
        current_x_counts <- colSums(current_x_counts)
        current_y_counts <- colSums(current_y_counts)
        out <- cbind(out, as.data.table(as.list(c(current_x_counts, current_y_counts))))
        out
    }
    collapsed_bin_counts <- realigned_block_counts[,collapse_bins(.SD, x_samples, y_samples),by=bin]

    ## get dp_matrix
    tumors.x <- paste0(tumor_samples,'.x')
    tumors.y <- paste0(tumor_samples,'.y')
    normal.x <- paste0(normal_sample,'.x')
    normal.y <- paste0(normal_sample,'.y')
    samples.x <- paste0(samples,'.x')
    samples.y <- paste0(samples,'.y')
    dp_matrix <- collapsed_bin_counts[,c(samples.x),with=F] + collapsed_bin_counts[,c(samples.y),with=F]
    colnames(dp_matrix) <- samples.y
    collapsed_bin_counts$dp=rowSums(dp_matrix)

    ## calculate BAF among bins with at least 10 reads in that sample
    baf_matrix <- collapsed_bin_counts[,c(samples.y),with=F] / dp_matrix
    baf_matrix[dp_matrix < min_block_reads] <- NA
    colnames(baf_matrix) <- samples
    res <- cbind(collapsed_bin_counts[,c('bin','bin_p','bin_q','dp'),with=F], baf_matrix) 
    resM <- melt(res, id.vars=c('bin','bin_p','bin_q','dp'))
    setnames(resM,c('variable','value'),c('sample','baf'))

    ## combine collapsed bin bafs with the count data
    d <- merge(cnt, resM, by=c('sample','bin'), all.x=T)
    d <- d[order(chr,bin_start,bin_end),]
    d$bin <- factor(d$bin, levels=unique(d$bin)) 
    dN <- d[sample==normal_sample,c('bin','count')]
    setnames(dN,'count','n_count')
    d <- merge(d, dN, by=c('bin'), all.x=T)

    ## calculate LogR
    .get_LogR <- function(d) {
        mid <- median(d$count[d$chr %in% c(1:22)],na.rm=T)
        d$LogR <- log2(d$count / mid)
        d
    }
    d <- d[,.get_LogR(.SD),by=sample]

    d[is.nan(baf),baf:=NA]
    d[!is.na(baf),count1:=count * baf]
    d[!is.na(baf),count2:=count * (1-baf)]
    d[is.na(baf), count1:=count]
    d[is.na(baf), count2:=0]
    d[is.na(baf),Ref:='.']
    d[is.na(baf),Alt:='.']
    d[!is.na(baf),Ref:='N']
    d[!is.na(baf),Alt:='N']
    d[,position:=round((bin_start + bin_end - 1)/2)]
    d <- d[,c('chr','position','Ref','Alt','sample','count1','count2','LogR'),with=F]
    d <- d[order(chr,position,sample),]

    message('Adding LogR and BAF values for each bin ...')
    d[,dp:=count1+count2]
    d[Ref=='N', BAF:=count1 / dp]
    d[BAF < 0 | BAF > 1, BAF:=NA]
    d[,bin:=paste0(chr,':',position)]
    d$bin <- as.integer(factor(d$bin, levels=unique(d$bin)))

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
    message('Done!')
    d
}
