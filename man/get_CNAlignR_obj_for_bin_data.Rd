% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_CNAlignR_obj_for_bin_data.R
\name{get_CNAlignR_obj_for_bin_data}
\alias{get_CNAlignR_obj_for_bin_data}
\title{get_CNAlignR_obj_for_bin_data}
\usage{
get_CNAlignR_obj_for_bin_data(
  qdnaseq_data,
  pileup_data,
  phased_bcf,
  sample_map,
  patient,
  sex,
  normal_sample,
  build,
  data_dir = ".",
  max_phaseable_distance = 20000,
  min_bin_reads_for_baf = 10,
  blacklisted_regions_file = NA,
  LogR_range_allowed = c(-3, 3),
  LogR_winsor_percentiles = c(NA, NA),
  LogR_smooth_bins = NA,
  normal_correction = F,
  multipcf_penalty = 70,
  multipcf_refine = F,
  multipcf_selectAlg = "exact",
  cleanup = T,
  seed = NA
)
}
\description{
Load and preprocess SCNA bin data from GLIMPSE2, QDNAseq and snp-pileup workflows.
}
