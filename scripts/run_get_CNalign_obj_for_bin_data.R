library(CNAlignR)

## define input arguments/parameters
phased_bcf <- 'glimpse/C149_G8_ligated_cleaned.bcf'
pileup_data <- 'C149_G8_ligated_cleaned.pileup.gz'
qdnaseq_data <- 'C149_100kbp_withXYM_hg38.rds'
sample_map <- 'C149_CNAlignR_sample_map.txt'
patient <- 'C149'
sex <- 'XX'
normal_sample <- 'G8'
build <- 'hg38'
data_dir <- '.'
seed=42
normal_correction=F
multipcf_penalty=70

## get CNalign data object for lowpass input data
obj <- get_CNAlignR_obj_for_bin_data(qdnaseq_data=qdnaseq_data,
                                    pileup_data=pileup_data,
                                    phased_bcf=phased_bcf,
                                    sample_map=sample_map,
                                    patient=patient,
                                    sex=sex,
                                    normal_sample=normal_sample,
                                    build=build,
                                    data_dir=data_dir,
                                    seed=seed,
                                    multipcf_penalty=multipcf_penalty,
                                    normal_correction=normal_correction)

## save the CNAlignR 'data object' to RDS file
saveRDS(obj, file=file.path(data_dir,paste0(patient,'_CNAlignR_obj.rds')))




