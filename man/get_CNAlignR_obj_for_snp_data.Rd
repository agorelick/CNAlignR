% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_CNAlignR_obj_for_snp_data.R
\name{get_CNAlignR_obj_for_snp_data}
\alias{get_CNAlignR_obj_for_snp_data}
\title{get_CNAlignR_obj_for_snp_data}
\usage{
get_CNAlignR_obj_for_snp_data(
  ascat_dir,
  sex,
  build,
  normal_sample,
  GCcontentfile = NULL,
  replictimingfile = NULL,
  multipcf_penalty = 300,
  multipcf_refine = F,
  multipcf_selectAlg = "exact",
  cleanup = T,
  seed = NA,
  output_dir = ".",
  obj_filename = "CNAlignR_obj.Rdata",
  Tumor_LogR_filename = "Tumor_LogR.txt",
  Tumor_BAF_filename = "Tumor_BAF.txt",
  Germline_LogR_filename = "Germline_LogR.txt",
  Germline_BAF_filename = "Germline_BAF.txt"
)
}
\description{
Load and preprocess standard SNP data (whole-exome, standard-depth whole genome). After running ASCAT on each tumor sample run individually (against the single normal sample), this script will load the resulting data for each sample and merge them, and run multi-sample segmentation via ASCAT's multipcf algorithm. It will then format the output for use with CNAlignR.
}
