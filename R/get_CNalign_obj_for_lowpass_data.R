##' get_CNalign_obj_for_lowpass_data
##'
##' Load and preprocess low-pass WGS data from GLIMPSE2, QDNAseq and snp-pileup workflows.
##'
##' @export
get_CNalign_obj_for_lowpass_data <- function(qdnaseq_data, pileup_data, phased_bcf, sample_map, patient, sex, normal_sample, build, data_dir='.', max_phaseable_distance=20000, min_bin_reads_for_baf=10, blacklisted_regions_file=NA, LogR_outlier_percentiles=c(NA,0.99),  multipcf_penalty=70, multipcf_refine=F, multipcf_selectAlg='exact', cleanup=F, seed=NA) {

    ## load and preprocess data from GLIMPSE2, QDNAseq and snp-pileup
    preprocessed_data <- preprocess_lowpass_data(qdnaseq_data=qdnaseq_data,
                                                 pileup_data=pileup_data,
                                                 phased_bcf=phased_bcf,
                                                 sample_map=sample_map,
                                                 normal_sample=normal_sample,
                                                 sex=sex,
                                                 build=build,
                                                 max_phaseable_distance=max_phaseable_distance,
                                                 min_bin_reads_for_baf=min_bin_reads_for_baf,
                                                 blacklisted_regions_file=blacklisted_regions_file,
                                                 LogR_outlier_percentiles=LogR_outlier_percentiles)

    ## define temporary file paths
    if(!dir.exists(data_dir)) dir.create(data_dir, recursive=T)
    Tumor_LogR_file <- file.path(data_dir,paste(patient,'Tumor_LogR.txt',sep='_'))
    Tumor_BAF_file <- file.path(data_dir,paste(patient,'Tumor_BAF.txt',sep='_'))
    Germline_LogR_file <- file.path(data_dir,paste(patient,'Germline_LogR.txt',sep='_'))
    Germline_BAF_file <- file.path(data_dir,paste(patient,'Germline_BAF.txt',sep='_'))

    ## create temporary files
    save_preprocessed_lowpass_data_for_ascat(preprocessed_data=preprocessed_data,
                                             normal_sample=normal_sample,
                                             sex=sex,
                                             Tumor_LogR_file=Tumor_LogR_file,
                                             Tumor_BAF_file=Tumor_BAF_file,
                                             Germline_LogR_file=Germline_LogR_file,
                                             Germline_BAF_file=Germline_BAF_file)

    ## get bin-level data
    message('Running preprocess_multipcf_data() ...')
    obj <- prep_data_for_multipcf(Tumor_LogR_file=Tumor_LogR_file,
                                  Tumor_BAF_file=Tumor_BAF_file,
                                  Germline_LogR_file=Germline_LogR_file,
                                  Germline_BAF_file=Germline_BAF_file,
                                  sex=sex,
                                  build=build)

    ## run multi-sample segmentation with default parameters
    message('Running run_ascat_multipcf() ...')
    obj <- run_ascat_multipcf(obj=obj,
                              build=build,
                              penalty=multipcf_penalty,
                              refine=multipcf_refine,
                              selectAlg=multipcf_selectAlg,
                              seed=seed)
            
    if(cleanup==T) {
        message('Cleaning up temporary files.')
        trash <- sapply(c(Tumor_LogR_file,Tumor_BAF_file,Germline_LogR_file,Germline_BAF_file), file.remove)
    }

    obj
}
 
