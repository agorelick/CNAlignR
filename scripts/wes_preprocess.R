library(ASCAT)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# prompts for command line arguments
parser$add_argument("--patient", type="character", help="Patient ID")
parser$add_argument("--tumor_name", type="character", help="Tumor sample name")
parser$add_argument("--tumor_bam", type="character", help="Path to tumor bam file")
parser$add_argument("--normal_name", type="character", help="Normal sample name")
parser$add_argument("--normal_bam", type="character", help="Path to normal bam file")
parser$add_argument("--sex", type="character", help="Patient sex, encoded as XX or XY")
parser$add_argument("--genome", type="character", help="Human reference genome version (hg19 or hg38, regardless of chr string in chromosome names)")
parser$add_argument("--target_bed", type="character", help="Path to .bed file with WES targeted regions")
parser$add_argument("--allelecounter_exe", type="character", help="Path to alleleCounter executable")
parser$add_argument("--alleles_prefix", type="character", help="Path+prefix for ASCAT alleles files (for each chromosome)")
parser$add_argument("--loci_prefix", type="character", help="Path+prefix for ASCAT loci files (for each chromosome)")
args <- parser$parse_args()

# return parsed arguments
cat("\nUser provided:\n")
cat("patient: ",args$patient,"\n")
cat("tumor_name: ",args$tumor_name,"\n")
cat("tumor_bam: ",args$tumor_bam,"\n")
cat("normal_name: ",args$normal_name,"\n")
cat("normal_bam: ",args$normal_bam,"\n")
cat("sex: ",args$sex,"\n")
cat("genome: ",args$genome,"\n")
cat("target_bed: ",args$target_bed,"\n")
cat("allelecounter_exe: ",args$allelecounter_exe,"\n")
cat("alleles_prefix: ",args$alleles_prefix,"\n")
cat("loci_prefix: ",args$loci_prefix,"\n")
cat("\n")

# output filenames
tumourLogR_file = paste(args$patient,args$tumor_name,args$normal_name,"Tumor_LogR.txt",sep='_')
tumourBAF_file = paste(args$patient,args$tumor_name,args$normal_name,"Tumor_BAF.txt",sep='_')
normalLogR_file = paste(args$patient,args$tumor_name,args$normal_name,"Germline_LogR.txt",sep='_')
normalBAF_file = paste(args$patient,args$tumor_name,args$normal_name,"Germline_BAF.txt",sep='_')

# run prepareHTS
ascat.prepareHTS(
                 tumourseqfile = args$tumor_bam,
                 tumourname = args$tumor_name,
                 normalseqfile = args$normal_bam,
                 normalname = paste0(args$tumor_name,'_',args$normal_name),
                 allelecounter_exe = args$allelecounter_exe,
                 alleles.prefix = args$alleles_prefix,
                 loci.prefix = args$loci_prefix,
                 gender = args$sex,
                 genomeVersion = args$genome,
                 nthreads = 1, # this MUST be 1 to prevent parallelization conflicts
                 tumourLogR_file = tumourLogR_file,
                 tumourBAF_file = tumourBAF_file,
                 normalLogR_file = normalLogR_file,
                 normalBAF_file = normalBAF_file,
                 BED_file = args$target_bed)


