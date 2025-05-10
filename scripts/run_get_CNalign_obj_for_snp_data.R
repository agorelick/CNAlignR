suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# essential parameters
parser$add_argument("--ascat_dir", type="character", help="Directory with output from alleleCounter")
parser$add_argument("--normal_sample", type="character", help="Normal sample name")

# output paths/filenames
parser$add_argument("--output_dir", type="character", help="output directory (where tmp files and data obj will go)", default='.')
parser$add_argument("--tumorlogr_file", type="character", help="Tumor LogR file", default='Tumor_LogR.txt')
parser$add_argument("--tumorbaf_file", type="character", help="Tumor BAF file", default='Tumor_BAF.txt')
parser$add_argument("--germlinelogr_file", type="character", help="Germline LogR file", default='Germline_LogR.txt')
parser$add_argument("--germlinebaf_file", type="character", help="Germline BAF file", default='Germline_BAF.txt')

# parse args
args <- parser$parse_args()

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
merged_tumorlogr_file <- file.path(output_dir,tumorlogr_file)
merged_tumorbaf_file <- file.path(output_dir,tumorbaf_file)
merged_germlinelogr_file <- file.path(output_dir,germlinelogr_file)
merged_germlinebaf_file <- file.path(output_dir,germlinebaf_file)
write.table(t_LogR, file = merged_tumorlogr_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_LogR, file = merged_germlinelogr_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(t_BAF, file = merged_tumorbaf_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(n_BAF, file = merged_germlinebaf_file, sep = "\t", quote = FALSE, col.names = NA)


