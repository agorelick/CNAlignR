library(CNalign)

## define input arguments/parameters
phased_bcf <- 'glimpse/SD28_NC1_ligated_cleaned.bcf'
pileup_data <- 'SD28_NC1_ligated_cleaned.pileup.gz'
qdnaseq_data <- 'SD28_100kbp_withXYM_hg38.rds'
sample_map <- 'SD_sample_map.txt'
patient <- 'SD28'
sex <- 'XX'
normal_sample <- 'NC1'
build <- 'hg38'
data_dir <- '.'
seed=42

## get CNalign data object for lowpass input data
obj <- get_CNalign_obj_for_bin_data(qdnaseq_data=qdnaseq_data,
                                    pileup_data=pileup_data,
                                    phased_bcf=phased_bcf,
                                    sample_map=sample_map,
                                    patient=patient,
                                    sex=sex,
                                    normal_sample=normal_sample,
                                    build=build,
                                    data_dir=data_dir,
                                    seed=seed)

## save the CNalign data object 'obj' to file (e.g., 'processed_data/C157_CNalign_obj.rds')
saveRDS(obj, file=file.path(data_dir,paste0(patient,'_CNalign_obj_no_normalization.rds')))




