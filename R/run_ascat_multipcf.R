##' run_ascat_multipcf
##' Run ascat.asmultipcf to do multi-sample segmentation, returns appended data objected.
##' @export
run_ascat_multipcf <- function(obj, penalty=70, seed=as.integer(Sys.time()), refine=T, selectAlg="exact", overwrite_mpcf=F, build='hg19') {
    chr_levels <- c(1:22,'X','Y')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Do multi-sample segmentation with ascat.asmultipcf
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(!'ascat.mpcf' %in% names(obj) | overwrite_mpcf==T) {
        message('Running ascat.asmultipcf() with penalty=',penalty,', seed=',seed,', selectAlg=',selectAlg,'.')
        ascat.mpcf = ascat.asmultipcf(obj$ascat.bc, out.dir=NA, seed=seed, penalty=penalty, refine=refine, selectAlg=selectAlg) 
        obj$ascat.mpcf <- ascat.mpcf
        obj$ascat.asmultipcf.params <- list(penalty=penalty, seed=seed, refine=refine, selectAlg=selectAlg)
    } else {
        message('obj already contains ascat.mpcf data, not overwriting it.')
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Extract segment data with start/end range and LogR/BAF values
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    message('Extracting segment data ...')

    ## extract the raw BAF and LogR from the ascat.mpcf object
    BAF_raw <- obj$ascat.mpcf$Tumor_BAF
    LogR_raw <- obj$ascat.mpcf$Tumor_LogR
    raw_common_rows <- intersect(rownames(BAF_raw), rownames(LogR_raw))
    BAF_raw <- BAF_raw[raw_common_rows,]
    LogR_raw <- LogR_raw[raw_common_rows,]

    ## extract the segmented LogR from the ascat.mpcf object
    LogR_segmented <- obj$ascat.mpcf$Tumor_LogR_segmented

    ## extract segmented BAF from ascat.mpcf object. then fill a matrix of the same size as the marker-level matrix with the non-NA BAFs
    BAF_segmented <- matrix(nrow=nrow(BAF_raw), ncol=ncol(BAF_raw))
    rownames(BAF_segmented) <- rownames(BAF_raw)
    colnames(BAF_segmented) <- colnames(BAF_raw)
    for(i in 1:ncol(BAF_segmented)) {
        vals <- obj$ascat.mpcf$Tumor_BAF_segmented[[i]] 
        BAF_segmented[rownames(vals),i] <- as.numeric(vals)
    }

    message('Extraploting from the non-NA BAFs to the rest of each chromosome ...')
    BAFs_avail <- BAF_segmented[!is.na(rowMeans(BAF_segmented)),]
    current_BAFs <- BAFs_avail[1,]
    #BAF_segmented <- extrapolate_BAFs(BAF_segmented, initial_BAFs) 

    for(i in 1:nrow(BAF_segmented)) {
        if(all(is.na(BAF_segmented[i,]))) {
            ## BAFs are currently NAs, use the last available value
            BAF_segmented[i,] <- current_BAFs
        } else {
            ## BAFs have non-NA values, use these and update the last available value
            current_BAFs <- BAF_segmented[i,] 
        }
    }

    message('Extracting SNPpos data from ascat.mpcf and add chr-arms to it ...')
    snp_pos <- obj$ascat.mpcf$SNPpos
    snp_pos <- cbind(VariantID=rownames(snp_pos), as.data.table(snp_pos))
    GD <- genome_data(build)
    arms <- GD$arms
    snp_pos[,pos2:=Position+1]
    snp_pos[,pos:=Position / 1e6]
    snp_pos[,pos2:=pos2 / 1e6]
    setkey(arms,'chr','arm_start','arm_end')
    setkey(snp_pos,'Chromosome','pos','pos2')
    snp_pos <- foverlaps(snp_pos, arms, type='within')
    snp_pos[,c('pos','pos2'):=NULL]
    snp_pos[,charm:=paste0(Chromosome,arm)]

    ## get number of germline copies (germline_copies) for each chromosome including sex chromosomes
    snp_pos[Chromosome %in% 1:22,germline_copies:=2]
    sex <- obj$ascat.loadData.params$gender
    if(sex=='XX') {
        snp_pos[Chromosome %in% 'X',germline_copies:=2]
        snp_pos[Chromosome %in% 'Y',germline_copies:=0]
    } else if(sex=='XY') {
        snp_pos[Chromosome %in% 'X',germline_copies:=1]
        snp_pos[Chromosome %in% 'Y',germline_copies:=1]   
    } else {
        snp_pos[Chromosome %in% c('X','Y'),germline_copies:=NA]
    }

    message('Finding unique segments based on any changes in BAF or LogR across any sample ...')
    current_segment <- 1
    current_BAFs <- BAF_segmented[1,]
    current_LogRs <- LogR_segmented[1,]
    segment <- rep(as.integer(NA), nrow(snp_pos))
    segment[1] <- current_segment
    current_charm <- snp_pos$charm[1]

    for(i in 2:nrow(snp_pos)) {
        new_BAFs <- BAF_segmented[i,]
        new_LogRs <- LogR_segmented[i,]
        new_charm <- snp_pos$charm[i]

        BAFs_avail <- which(!is.na(new_BAFs) & !is.na(current_BAFs))        
        BAF_changed <- any(new_BAFs[BAFs_avail] != current_BAFs[BAFs_avail])
        LogRs_avail <- which(!is.na(new_LogRs) & !is.na(current_LogRs))        
        LogR_changed <- any(new_LogRs[LogRs_avail] != current_LogRs[LogRs_avail])
        arm_changed <- new_charm!=current_charm

        if(BAF_changed | LogR_changed | arm_changed) { 
            current_segment <- current_segment + 1
            current_BAFs <- new_BAFs
            current_LogRs <- new_LogRs 
            current_charm <- new_charm
        } 
        segment[i] <- current_segment
    }
    snp_pos$segment <- segment
    
    message('Collapsing each segment to start/end ranges ...')
    get_segments <- function(snp_pos) {
        n_markers <- nrow(snp_pos)
        start <- min(snp_pos$Position)
        end <- max(snp_pos$Position)
        list(seg_start=start, seg_end=end, n_markers=n_markers)
    }
    segs <- snp_pos[,get_segments(.SD), by=c('Chromosome','arm','segment','germline_copies')]
    segs <- segs[order(segment),]
    segs$Chromosome <- factor(segs$Chromosome, levels=chr_levels)
    segs[,seg_length:=seg_end - seg_start + 1]
    if(any(duplicated(segs$segment))) stop('something went wrong.')
    obj$segments <- segs

    ## print out summary stats about segs here...
    message('N copy number segments: ',nrow(obj$segment))
   
    tbl <- quantile(obj$segments$seg_length,seq(0,1,by=0.1))
    tbl <- data.table(percentile=names(tbl), length_mb=round(as.numeric(tbl)/1e6,3)) 
    tbl <- as.data.frame(tbl[order(length_mb,decreasing=T),])
    message('Segment length deciles: ')
    print(tbl)


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Annotate the marker-level data with segment info, keep all markers
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    message('Annotating marker-level data with segment info ...')

    ## now annotate the complete snp_pos data with segment start/end info
    snp_pos[,Position2:=Position+1]
    setkey(snp_pos,'Chromosome','Position','Position2')
    setkey(segs,'Chromosome','seg_start','seg_end')
    snp_pos <- foverlaps(snp_pos, segs[,c('Chromosome','seg_start','seg_end','seg_length','n_markers'),with=F], type='any')
    snp_pos[,c('seg_start','seg_end','Position2','arm_start','arm_end'):=NULL]
        
    ## get long data.table with SNP-level RAW LogR and BAFs
    BAF_dat <- cbind(VariantID=rownames(BAF_raw), as.data.table(BAF_raw))
    BAF_dat <- merge(snp_pos, BAF_dat, by='VariantID', all.x=T)
    BAF_long <- data.table::melt(BAF_dat, id.vars=c('VariantID','Chromosome','segment','seg_length','n_markers','arm','Position','charm','germline_copies'))
    setnames(BAF_long,c('variable','value'),c('sample','BAF'))
    BAF_long$Chromosome <- factor(BAF_long$Chromosome, levels=chr_levels)
    BAF_long <- BAF_long[order(Chromosome,Position,sample),]   
    LogR_dat <- cbind(VariantID=rownames(LogR_raw), as.data.table(LogR_raw))
    LogR_dat <- merge(snp_pos, LogR_dat, by='VariantID', all.x=T)
    LogR_long <- data.table::melt(LogR_dat, id.vars=c('VariantID','Chromosome','segment','seg_length','n_markers','arm','Position','charm','germline_copies'))
    setnames(LogR_long,c('variable','value'),c('sample','LogR'))
    LogR_long$Chromosome <- factor(LogR_long$Chromosome, levels=chr_levels)
    LogR_long <- LogR_long[order(Chromosome,Position,sample),]   
    dat_long_raw <- cbind(LogR_long, BAF=BAF_long$BAF)  ## need to add error checking or do this less-efficiently by merging

    ## add segment-level LogR and BAFs
    LogR_dat <- cbind(VariantID=rownames(LogR_segmented), as.data.table(LogR_segmented))
    LogR_dat <- merge(snp_pos[,c('VariantID','segment'),with=F], LogR_dat, by='VariantID', all.y=T)
    LogR_dat <- LogR_dat[!duplicated(segment),]
    LogR_dat <- LogR_dat[order(segment),]
    LogR_dat[,VariantID:=NULL]
    LogR_long <- data.table::melt(LogR_dat, id.vars='segment')
    setnames(LogR_long,c('variable','value'),c('sample','LogR_segmented'))
    dat_long_raw <- merge(dat_long_raw, LogR_long, by=c('segment','sample'), all.x=T)    
    BAF_dat <- cbind(VariantID=rownames(BAF_segmented), as.data.table(BAF_segmented))
    BAF_dat <- merge(snp_pos[,c('VariantID','segment'),with=F], BAF_dat, by='VariantID', all.y=T)
    BAF_dat <- BAF_dat[!duplicated(segment),]
    BAF_dat <- BAF_dat[order(segment),]
    BAF_dat[,VariantID:=NULL]
    BAF_long <- data.table::melt(BAF_dat, id.vars='segment')
    setnames(BAF_long,c('variable','value'),c('sample','BAF_segmented'))
    dat_long_raw <- merge(dat_long_raw, BAF_long, by=c('segment','sample'), all.x=T)    

    ## append the global start pos (in Mb) to dat_long_raw
    chr_global_start <- GD$chr[,c('chr','global_start'),with=F]
    setnames(chr_global_start,'global_start','global_start_mb')
    dat_long_raw <- merge(dat_long_raw, chr_global_start, by.x='Chromosome', by.y='chr', all.x=T)
    dat_long_raw$Chromosome <- factor(dat_long_raw$Chromosome, levels=chr_levels)
    dat_long_raw <- dat_long_raw[order(Chromosome, Position, sample),]
    dat_long_raw[,global_pos_mb:=Position/1e6 + global_start_mb]
    dat_long_raw[,global_start_mb:=NULL]

    ## filter out BAF values from positions that are not likely het-SNPs based on the normal, because these will mess up the BAF ECDF stats later on
    tmp <- obj$ascat.bc$Germline_BAF
    tmp <- cbind(VariantID=rownames(tmp), as.data.table(tmp))
    names(tmp)[2] <- 'N'
    het_snps <- tmp[N >= 0.1 & N <= 0.9,(VariantID)] # this includes all chromosomes
    dat_long_raw[!VariantID %in% het_snps,BAF:=NA]
    obj$marker_level_annotated <- split(dat_long_raw, by='sample')


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Prepare segment-level BAF/LogRs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    message('Preparing segment-level data ...')

    LogR_segs <- merge(segs, LogR_dat, by='segment', all.x=T)
    LogR_segs <- merge(LogR_segs, chr_global_start, by.x='Chromosome', by.y='chr', all.x=T)
    LogR_segs[,global_seg_start_mb:=seg_start/1e6 + global_start_mb]
    LogR_segs[,global_seg_end_mb:=seg_end/1e6 + global_start_mb]
    LogR_segs[,global_start_mb:=NULL]
    LogR_segs <- data.table::melt(LogR_segs, id.vars=c('Chromosome','segment','arm','germline_copies','seg_start',
                                                       'seg_end','n_markers','seg_length','global_seg_start_mb','global_seg_end_mb'))
    setnames(LogR_segs,c('variable','value'),c('sample','LogR_segmented'))

    BAF_segs <- merge(segs, BAF_dat, by='segment', all.x=T)
    BAF_segs <- merge(BAF_segs, chr_global_start, by.x='Chromosome', by.y='chr', all.x=T)
    BAF_segs[,global_seg_start_mb:=seg_start/1e6 + global_start_mb]
    BAF_segs[,global_seg_end_mb:=seg_end/1e6 + global_start_mb]
    BAF_segs[,global_start_mb:=NULL]
    BAF_segs <- data.table::melt(BAF_segs, id.vars=c('Chromosome','segment','arm','germline_copies','seg_start',
                                                       'seg_end','n_markers','seg_length','global_seg_start_mb','global_seg_end_mb'))
    setnames(BAF_segs,c('variable','value'),c('sample','BAF_segmented'))

    ## add segment-level data to obj
    segment_level <- cbind(LogR_segs, BAF_segmented=BAF_segs$BAF_segmented)
    obj$segment_level <- split(segment_level, by='sample')


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # generate data formatted for CNalign
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ## long table of sample,segment-level BAF, LogR, GC for CNalign input
    CNalign_dat <- as.data.frame(segment_level[,c('sample','segment','LogR_segmented','BAF_segmented','germline_copies'),with=F])
    names(CNalign_dat) <- c('sample','segment','LogR','BAF','GC')
    obj$CNalign_dat <- CNalign_dat

 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Get ECDFs for LogR and BAFs for each sample
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    obj$fit_segments <- sort(unique(obj$segments[obj$segments$germline_copies==2 & obj$segments$seg_length >= 5e6,(segment)]))
    get_ecdf_per_sample <- function(ml) { 
        ml <- ml[!is.na(segment)]
        ml <- ml[segment %in% obj$fit_segments,]
        ml[,diff_from_BAF:=(BAF - BAF_segmented)^2]  
        ml[,diff_from_LogR:=(LogR - LogR_segmented)^2]
        LogR_CDF <- ecdf(ml$diff_from_LogR[!is.na(ml$diff_from_LogR)])
        BAF_CDF <- ecdf(ml$diff_from_BAF[!is.na(ml$diff_from_BAF)])
        list(LogR_CDF=LogR_CDF, BAF_CDF=BAF_CDF)
    }
    ECDF_fits <- lapply(obj$marker_level_annotated, get_ecdf_per_sample)
    obj$ECDF_fits <- ECDF_fits

    # return the result
    message('Done!')
    obj
}
