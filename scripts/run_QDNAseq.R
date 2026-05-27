### created segmented copy number profiles to include X and Y chromosome

library(data.table)
library(QDNAseq.hg38)
library(Biobase)

patient <- 'SD28'

bams <- c('original_bams/SD28_NC1_T1.bam', 'original_bams/SD28_NC2_T1.bam', 'original_bams/SD28_NC3_T1.bam', 'original_bams/SD28_NC4_T1.bam', 'original_bams/SD28_NC5_T1.bam', 'original_bams/SD28_NC6_T1.bam', 'original_bams/SD28_TC1-1-Rep1_T1.bam', 'original_bams/SD28_TC1-2_T1.bam', 'original_bams/SD28_TC1-3_T1.bam', 'original_bams/SD28_TC1-4_T1.bam', 'original_bams/SD28_TC2-1_T1.bam', 'original_bams/SD28_TC2-2_T1.bam')

get_corrected_bin_counts_for_sample <- function(bam, binsize) {
    message('\nRunning QDNAseq on ',bam)   

    ## file with bin-level annotations. Bins are 1Mb, but include chrM. Variable bin size is controlled by the 'bases' annotation (% of 1Mb with non-N base, which is 16,569/1e6 for chrM)
    suffix <- tail(strsplit(gsub('[.]bam','',bam),'[/]')[[1]],1)

    message('Loading bin information ...')
    ref_bins <- getBinAnnotations(binSize=binsize, genome='hg38')
    bin_dat <- as.data.frame(ref_bins@data)
    bin_dat <- cbind(feature=rownames(bin_dat), as.data.table(bin_dat))

    ## generating bin annotations of size 1MB
    message('binReadCounts ...')
    readCounts <- binReadCounts(ref_bins, bamfiles=bam) # for each bam file

    ## Make sure sex chromosomes are excluded here, so that we do not include them in calculating the loess correction
    message('applyFilters ...')
    readCountsFiltered <- applyFilters(readCounts, residual=T, blacklist=T, mappability=25, chromosomes=c('X','Y','MT')) ## filters out X,Y,MT before estimating correction factors

    ## do loess correction
    message('estimateCorrection ...')
    readCountsFiltered <- estimateCorrection(readCountsFiltered)

    ## now we add X,Y back in and apply the correction that we obtained previously 
    message('applyFilters (adding X, Y, MT) ...')
    readCountsFiltered <- applyFilters(readCountsFiltered, residual=T, blacklist=T, mappability=25, chromosomes=c()) 
    message('correctBins ...')
    correctedBins <- correctBins(readCountsFiltered)  

    ## continue with standard QDNAseq pipeline
    correctedBinsNormalized <- normalizeBins(correctedBins)

    ## export the bin data to a TSV so that we can read it back and do some formatting, then re-save it
    tmpfilenormalized <- paste0('~/scratch/',suffix,'_',binsize,'kbp_withXYM_hg38_correctednormalized.txt')
    exportBins(correctedBinsNormalized, file=tmpfilenormalized, format='tsv', logTransform=F, type='copynumber')
    bins_correctednormalized <- fread(tmpfilenormalized)
    names(bins_correctednormalized)[ncol(bins_correctednormalized)] <- 'count'
    out <- bins_correctednormalized
    out$chromosome <- factor(out$chromosome, levels=c(1:22,'X','Y','MT'))
    out <- out[!is.na(chromosome) & !is.na(start),]
    out <- out[order(chromosome, start, end),]
    out$file <- bam
    out$binsize_kbp <- binsize
    out$genome <- 'hg38'

    out
}

l <- lapply(bams, get_corrected_bin_counts_for_sample, binsize=100)
out <- rbindlist(l)
saveRDS(out,file=paste0(patient,'_100kbp_withXYM_hg38.rds'))




