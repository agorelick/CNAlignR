library(CNalign)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# essential parameters
parser$add_argument("--ascat_dir", type="character", help="Directory with output from alleleCounter")
parser$add_argument("--patient", type="character", help="Patient ID")
parser$add_argument("--sex", type="character", help="Patient sex, encoded as XX or XY")
parser$add_argument("--genome", type="character", help="Human reference genome version (hg19 or hg38, regardless of chr string in chromosome names)")
parser$add_argument("--normal_name", type="character", help="Normal sample name")
parser$add_argument("--gc_file", type="character", help="GCcontentfile", default=NULL)
parser$add_argument("--rt_file", type="character", help="replictimingfile", default=NULL)

# output paths/filenames
parser$add_argument("--output_dir", type="character", help="output directory (where tmp files and data obj will go)", default='.')
parser$add_argument("--obj_file", type="character", help="CNalign data object filename", default='CNalign_obj.Rdata')
parser$add_argument("--tumorlogr_file", type="character", help="Tumor LogR file", default='Tumor_LogR.txt')
parser$add_argument("--tumorbaf_file", type="character", help="Tumor BAF file", default='Tumor_BAF.txt')
parser$add_argument("--germlinelogr_file", type="character", help="Germline LogR file", default='Germline_LogR.txt')
parser$add_argument("--germlinebaf_file", type="character", help="Germline BAF file", default='Germline_BAF.txt')

args <- parser$parse_args()

# run get_CNalign_obj_for_snp_data
CNalign::get_CNalign_obj_for_snp_data(ascat_dir=args$ascat_dir, 
                                      sex=args$sex, 
                                      build=args$genome, 
                                      normal_sample=args$normal_name,
                                      GCcontentfile=args$gc_file, 
                                      replictimingfile=args$rt_file, 
                                      output_dir=args$output_dir, 
                                      obj_filename=args$obj_file,
                                      Tumor_LogR_filename=args$tumorlogr_file,
                                      Tumor_BAF_filename=args$tumorbaf_file,
                                      Germline_LogR_filename=args$germlinelogr_file,
                                      Germline_BAF_filename=args$germlinebaf_file)


