##' get_CNalign_obj_for_ascat_data
##'
##' Load and preprocess ASCAT data (whole-exome, standard-depth whole genome).
##'
##' @export
get_CNalign_obj_for_ascat_data <- function(sample_map, patient, sex, normal_sample, build, data_dir='.', max_phaseable_distance=20000, min_bin_reads_for_baf=10, blacklisted_regions_file=NA, LogR_range_allowed=c(-3.0,3.0), LogR_winsor_percentiles=c(NA,NA), LogR_smooth_bins=NA, multipcf_penalty=70, multipcf_refine=F, multipcf_selectAlg='exact', cleanup=T, seed=NA) {

# After running ASCAT on each tumor sample run individually (against the single normal sample)
# this script will load the resulting data for each sample and merge them, and run multi-sample segmentation
# from ASCAT's multipcf algorithm.
# It will then format the output for use with CNalign.
#
# This script should be called in a slurm job for each patient with input argumentts for the respective patient

library(data.table)
library(ASCAT)
library(lpASCN)
source('/n/data1/hms/genetics/naxerova/lab/alex/cholangio_RA_bams_Broad/scripts/run_multisample_segmentation_functions.R') # source required functions (eventually put this in lpASCN package)

start_time <- Sys.time()
message('Starting at ', as.character(start_time))


## get input arguments
args = as.character(commandArgs(trailingOnly=TRUE))

ascat_dir <- args[1] # directory with ASCAT Tumor/Germline LogR/BAF for each sample
output_dir <- args[2] # directory for the outut destination
map_file <- args[3] # path to sample map file
patient <- args[4] # patient ID (to subset the map file for the right samples)
sex <- args[5] # patient sex coded as "XX" or "XY"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run for each CCA patient on O2 (default data processing workflow)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

map <- fread(map_file)
map <- map[Patient_ID %in% patient]
normal_sample <- map[Sample_Type=='Normal',(Sample_ID)]
tumor_samples <- map[Sample_Type=='Tumor',(Sample_ID)]
map[Sample_Type=='Tumor',prefix:=file.path(ascat_dir,paste(sample_name_prefix,Sample_ID,normal_sample,sep='_'))]
prefices <- map[!is.na(prefix),(prefix)]

get_data <- function(prefix, suffix) {
    file <- paste0(prefix,suffix)
    if(file.exists(file)) {
        message('Loading data from ',file)
        out <- fread(file)
        sample <- names(out)[4]
        out$sample <- sample
        names(out)[4] <- 'value'
        out$suffix <- suffix
        out
    } else {
        message(file,' does not exist!')
        NULL
    }
}

## get LogR data 
l_Tumor_LogR <- lapply(prefices, get_data, suffix='_Tumor_LogR.txt')
Tumor_LogR <- rbindlist(l_Tumor_LogR)
Tumor_LogR$Chromosome <- factor(Tumor_LogR$Chromosome, levels=c(1:22,'X','Y'))
l_Germline_LogR <- lapply(prefices, get_data, suffix='_Germline_LogR.txt')
Germline_LogR <- rbindlist(l_Germline_LogR)
Germline_LogR <- Germline_LogR[!duplicated(V1),]
Germline_LogR[,sample:=normal_sample]
Germline_LogR$Chromosome <- factor(Germline_LogR$Chromosome, levels=c(1:22,'X','Y'))
x_LogR <- rbind(Tumor_LogR, Germline_LogR)
x_LogR <- data.table::dcast(V1 + Chromosome + Position ~ sample, value.var='value', data=x_LogR)
x_LogR <- x_LogR[order(Chromosome, Position),]
x_LogR <- as.data.frame(x_LogR)
rownames(x_LogR) <- x_LogR$V1
x_LogR$V1 <- NULL

## get BAF data 
l_Tumor_BAF <- lapply(prefices, get_data, suffix='_Tumor_BAF.txt')
Tumor_BAF <- rbindlist(l_Tumor_BAF)
Tumor_BAF$Chromosome <- factor(Tumor_BAF$Chromosome, levels=c(1:22,'X','Y'))
l_Germline_BAF <- lapply(prefices, get_data, suffix='_Germline_BAF.txt')
Germline_BAF <- rbindlist(l_Germline_BAF)
Germline_BAF <- Germline_BAF[!duplicated(V1),]
Germline_BAF[,sample:=normal_sample]
Germline_BAF$Chromosome <- factor(Germline_BAF$Chromosome, levels=c(1:22,'X','Y'))
x_BAF <- rbind(Tumor_BAF, Germline_BAF)
x_BAF <- data.table::dcast(V1 + Chromosome + Position ~ sample, value.var='value', data=x_BAF)
x_BAF <- x_BAF[order(Chromosome, Position),]
x_BAF <- as.data.frame(x_BAF)
rownames(x_BAF) <- x_BAF$V1
x_BAF$V1 <- NULL

## subset for positions in both logR and BAF data
valid_rows <- intersect(rownames(x_BAF), rownames(x_LogR))
x_LogR <- x_LogR[valid_rows,]
x_BAF <- x_BAF[valid_rows,]
t_LogR <- x_LogR[,c('Chromosome','Position',tumor_samples)]
t_BAF <- x_BAF[,c('Chromosome','Position',tumor_samples)]
n_LogR <- x_LogR[,c('Chromosome','Position',normal_sample)]
n_BAF <- x_BAF[,c('Chromosome','Position',normal_sample)]

#if(!dir.exists(output_dir)) dir.create(output_dir)
tumor_LogR_file <- file.path(output_dir,'Tumor_LogR.txt')
tumor_BAF_file <- file.path(output_dir,'Tumor_BAF.txt')
germline_LogR_file <- file.path(output_dir,'Germline_LogR.txt')
germline_BAF_file <- file.path(output_dir,'Germline_BAF.txt')
write.table(t_LogR, file = tumor_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_LogR, file = germline_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(t_BAF, file = tumor_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_BAF, file = germline_BAF_file, sep = "\t", quote = FALSE, col.names = NA)

## get bin-level data
message('Running preprocess_multipcf_data() ...')
obj <- preprocess_multipcf_data(Tumor_LogR_file = tumor_LogR_file,
                          Tumor_BAF_file = tumor_BAF_file,
                          Germline_LogR_file = germline_LogR_file,
                          Germline_BAF_file = germline_BAF_file,
                          sex = sex,
                          genomeVersion = 'hg19',
                          GCcontentfile = '/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/GC_G1000_hg19.txt',
                          replictimingfile='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/RT_G1000_hg19.txt')

message('done! saving result before trying segmentation ...')
saveRDS(obj, file=file.path(output_dir,paste0(patient,'_CNalign_tmpobj.rds')))


## add segment-level data
message('Running run_ascat_multipcf() ...')
obj <- run_ascat_multipcf(obj, seed=42)
saveRDS(obj, file=file.path(output_dir,paste0(patient,'_CNalign_obj.rds')))
message('Done!')

end_time <- Sys.time()
message('Finished at ', as.character(end_time))


