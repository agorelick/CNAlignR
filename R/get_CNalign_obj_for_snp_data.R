##' get_CNalign_obj_for_snp_data
##'
##' Load and preprocess standard SNP data (whole-exome, standard-depth whole genome). After running ASCAT on each tumor sample run individually (against the single normal sample), this script will load the resulting data for each sample and merge them, and run multi-sample segmentation via ASCAT's multipcf algorithm. It will then format the output for use with CNalign.
##'
##' @export
get_CNalign_obj_for_snp_data <- function(sample_map, patient, sex, build, ascat_dir='.', data_dir='.', GCcontentfile=NULL, replictimingfile=NULL, multipcf_penalty=70, multipcf_refine=F, multipcf_selectAlg='exact', cleanup=T, seed=NA) {

    map <- fread(sample_map)
    map <- map[patient_id %in% patient]
    normal_sample <- map[sample_type=='Normal',(lesion_id)]
    map <- map[sample_type=='Tumor']
    tumor_samples <- map[sample_type=='Tumor',(lesion_id)]
    prefices <- file.path(ascat_dir,map[!is.na(ascat_prefix),(ascat_prefix)])

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
    tumor_samples <- tumor_samples[tumor_samples %in% names(x_LogR)]
    valid_rows <- intersect(rownames(x_BAF), rownames(x_LogR))
    x_LogR <- x_LogR[valid_rows,]
    x_BAF <- x_BAF[valid_rows,]
    t_LogR <- x_LogR[,c('Chromosome','Position',tumor_samples)]
    t_BAF <- x_BAF[,c('Chromosome','Position',tumor_samples)]
    n_LogR <- x_LogR[,c('Chromosome','Position',normal_sample)]
    n_BAF <- x_BAF[,c('Chromosome','Position',normal_sample)]

    if(!dir.exists(data_dir)) dir.create(data_dir)
    Tumor_LogR_file <- file.path(data_dir,'Tumor_LogR.txt')
    Tumor_BAF_file <- file.path(data_dir,'Tumor_BAF.txt')
    Germline_LogR_file <- file.path(data_dir,'Germline_LogR.txt')
    Germline_BAF_file <- file.path(data_dir,'Germline_BAF.txt')
    write.table(t_LogR, file = Tumor_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_LogR, file = Germline_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(t_BAF, file = Tumor_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_BAF, file = Germline_BAF_file, sep = "\t", quote = FALSE, col.names = NA)

    ## get snp-level data
    message('Running preprocess_multipcf_data() ...')
    obj <- prep_data_for_multipcf(Tumor_LogR_file=Tumor_LogR_file,
                                  Tumor_BAF_file=Tumor_BAF_file,
                                  Germline_LogR_file=Germline_LogR_file,
                                  Germline_BAF_file=Germline_BAF_file,
                                  sex=sex,
                                  build=build,
                                  GCcontentfile=GCcontentfile,
                                  replictimingfile=replictimingfile)

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

    main_params <- list(
                   sample_map=sample_map,
                   patient=patient,
                   sex=sex,
                   normal_sample=normal_sample,
                   build=build,
                   ascat_dir=ascat_dir,
                   data_dir=data_dir,
                   GCcontentfile=GCcontentfile,
                   replictimingfile=replictimingfile,
                   multipcf_penalty=multipcf_penalty,
                   multipcf_refine=multipcf_refine,
                   multipcf_selectAlg=multipcf_selectAlg,
                   cleanup=cleanup,
                   seed=seed) 

    obj$main_params <- main_params
    obj
}
 

