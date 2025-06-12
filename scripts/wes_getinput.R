message('Starting script to create CNAlignR data object.')

message('Loading required packages ...')
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("CNAlignR"))

# create parser object
parser <- ArgumentParser()

# input parameters (required)
parser$add_argument("--patient", type="character", help="Patient ID")
parser$add_argument("--normal_sample", type="character", help="Normal sample name")
parser$add_argument("--sex", type="character", help="Patient sex, encoded XX or XY")
parser$add_argument("--build", type="character", help="Human reference genome version, encoded hg19 or hg38 (regardless of chr prefix in contig names)")
parser$add_argument("--GCcontentfile", type="character", help="GCcontentfile")
parser$add_argument("--replictimingfile", type="character", help="replictimingfile")

# multipcf arguments
parser$add_argument("--penalty", type="numeric", help="Penalty for multipcf, should be large (~300) for purity/ploidy fitting", default=300)
parser$add_argument("--penalty_hisens", type="numeric", help="Penalty for high-sensitivity multipcf, should be small (~25) for CCF calculations", default=25)
parser$add_argument("--seed", type="integer", help="seed for multipcf (-1 = no seed)", default=-1)
parser$add_argument("--selectAlg", type="character", help="selectAlg for multipcf (exact|fast)", default='exact')
parser$add_argument("--refine", type="logical", help="should we refine segments after multipcf?", default=TRUE)

# input parameters (optional)
parser$add_argument("--ascat_dir", type="character", help="Directory with output from alleleCounter", default='.')
parser$add_argument("--output_dir", type="character", help="output directory (where tmp files and data obj will go)", default='.')
parser$add_argument("--tumorlogr_file", type="character", help="Tumor LogR file", default='Tumor_LogR.txt')
parser$add_argument("--tumorbaf_file", type="character", help="Tumor BAF file", default='Tumor_BAF.txt')
parser$add_argument("--germlinelogr_file", type="character", help="Germline LogR file", default='Germline_LogR.txt')
parser$add_argument("--germlinebaf_file", type="character", help="Germline BAF file", default='Germline_BAF.txt')
parser$add_argument("--obj_file", type="character", help="CNalign data object file", default='CNalign_obj.rds')
parser$add_argument("--overwrite", type="logical", help="Overwrite existing data? (TRUE|FALSE)", default=FALSE)


# parse args
args <- parser$parse_args()
patient = args$patient
normal_sample = args$normal_sample
sex = args$sex
build = args$build
GCcontentfile = args$GCcontentfile
replictimingfile = args$replictimingfile
ascat_dir = args$ascat_dir
output_dir = args$output_dir
tumorlogr_file = args$tumorlogr_file
tumorbaf_file = args$tumorbaf_file
germlinelogr_file = args$germlinelogr_file
germlinebaf_file = args$germlinebaf_file
obj_file = args$obj_file
overwrite = args$overwrite


# dynamically get input file names
ascat_files <- dir(ascat_dir,full.names=T)
Tumor_LogR_files <- grep(paste0(normal_sample,'_Tumor_LogR.txt'), ascat_files, value=T)
Tumor_BAF_files <- grep(paste0(normal_sample,'_Tumor_BAF.txt'), ascat_files, value=T)
Germline_LogR_files <- grep(paste0(normal_sample,'_Germline_LogR.txt'), ascat_files, value=T)
Germline_BAF_files <- grep(paste0(normal_sample,'_Germline_BAF.txt'), ascat_files, value=T)

# helper to load data file
get_data <- function(file) {
    if(file.exists(file)) {
        message('Loading data from ',file)
        out <- fread(file)
        sample <- names(out)[4]
        out$sample <- sample
        names(out)[4] <- 'value'
        out
    } else {
        message(file,' does not exist!')
        NULL
    }
}

# get LogR data 
l_Tumor_LogR <- lapply(Tumor_LogR_files, get_data)
Tumor_LogR <- rbindlist(l_Tumor_LogR)
Tumor_LogR$Chromosome <- factor(Tumor_LogR$Chromosome, levels=c(1:22,'X','Y'))
l_Germline_LogR <- lapply(Germline_LogR_files, get_data)
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

# get BAF data 
l_Tumor_BAF <- lapply(Tumor_BAF_files, get_data)
Tumor_BAF <- rbindlist(l_Tumor_BAF)
Tumor_BAF$Chromosome <- factor(Tumor_BAF$Chromosome, levels=c(1:22,'X','Y'))
l_Germline_BAF <- lapply(Germline_BAF_files, get_data)
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

# subset for positions in both logR and BAF data
all_samples <- names(x_LogR)[!names(x_LogR) %in% c('V1','Chromosome','Position')]
tumor_samples <- all_samples[!all_samples %in% normal_sample]
valid_rows <- intersect(rownames(x_BAF), rownames(x_LogR))
x_LogR <- x_LogR[valid_rows,]
x_BAF <- x_BAF[valid_rows,]
t_LogR <- x_LogR[,c('Chromosome','Position',tumor_samples)]
t_BAF <- x_BAF[,c('Chromosome','Position',tumor_samples)]
n_LogR <- x_LogR[,c('Chromosome','Position',normal_sample)]
n_BAF <- x_BAF[,c('Chromosome','Position',normal_sample)]

# save the merged data
if(!dir.exists(output_dir)) dir.create(output_dir)
merged_tumorlogr_file <- file.path(output_dir,paste0(patient,'_',tumorlogr_file))
merged_tumorbaf_file <- file.path(output_dir,paste0(patient,'_',tumorbaf_file))
merged_germlinelogr_file <- file.path(output_dir,paste0(patient,'_',germlinelogr_file))
merged_germlinebaf_file <- file.path(output_dir,paste0(patient,'_',germlinebaf_file))
write.table(t_LogR, file = merged_tumorlogr_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_LogR, file = merged_germlinelogr_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(t_BAF, file = merged_tumorbaf_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_BAF, file = merged_germlinebaf_file, sep = "\t", quote = FALSE, col.names = NA)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prep data for multipcf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output_path <- file.path(output_dir, obj_file)
if(!file.exists(output_path) | overwrite==T) { 
    ## get snp-level data
    message('Running prep_data_for_multipcf() ...')
    obj <- prep_data_for_multipcf(Tumor_LogR_file=merged_tumorlogr_file,
                                  Tumor_BAF_file=merged_tumorbaf_file,
                                  Germline_LogR_file=merged_germlinelogr_file,
                                  Germline_BAF_file=merged_germlinebaf_file,
                                  sex=sex,
                                  build=build,
                                  GCcontentfile=GCcontentfile,
                                  replictimingfile=replictimingfile)

    ## add a list with the parameters to the output object so that we have a record of what values were used
    main_params <- list(patient = args$patient,
                        normal_sample = args$normal_sample,
                        sex = args$sex,
                        build = args$build,
                        GCcontentfile = args$GCcontentfile,
                        replictimingfile = args$replictimingfile,
                        ascat_dir = args$ascat_dir,
                        output_dir = args$output_dir,
                        tumorlogr_file = args$tumorlogr_file,
                        tumorbaf_file = args$tumorbaf_file,
                        germlinelogr_file = args$germlinelogr_file,
                        germlinebaf_file = args$germlinebaf_file,
                        obj_file = args$obj_file)
    obj$main_params <- main_params

    # save object 
    message('Saving CNalign data object to: ', output_path)
    saveRDS(obj, file=output_path)
} else {
    message(output_path,' already exists. Not overwriting it.')
    obj <- readRDS(output_path)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run multipcf for purity/ploidy fitting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output2_path <- gsub('[.]rds','_mpcf.rds',output_path)
if(!file.exists(output2_path) | overwrite==T) { 
    # do multipcf with large penalty (better for purity/ploidy fitting)
    myseed = ifelse(args$seed!=-1,args$seed,as.numeric(NA))
    obj2 <- run_ascat_multipcf(obj, build=args$build, penalty=args$penalty, seed=myseed, selectAlg=args$selectAlg, refine=args$refine)
    obj2$main_params <- list(build=args$build, penalty=args$penalty, seed=myseed, selectAlg=args$selectAlg, refine=args$refine)
    saveRDS(obj2, file=output2_path)
} else {
    message(output2_path,' already exists. Not overwriting it.')
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# do hisens multipcf with small penalty
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output3_path <- gsub('[.]rds','_mpcf_hisens.rds',output_path)
if(!file.exists(output3_path) | overwrite==T) { 
    myseed = ifelse(args$seed!=-1,args$seed,as.numeric(NA))
    obj3 <- run_ascat_multipcf(obj, build=args$build, penalty=args$penalty_hisens, seed=myseed, selectAlg=args$selectAlg, refine=args$refine)
    obj3$main_params <- list(build=args$build, penalty_hisens=args$penalty_hisens, seed=myseed, selectAlg=args$selectAlg, refine=args$refine)
    saveRDS(obj3, file=output3_path)
} else {
    message(output3_path,' already exists. Not overwriting it.')
}




