##' get_CNAlignR_obj_for_snp_data
##'
##' Load and preprocess standard SNP data (whole-exome, standard-depth whole genome). After running ASCAT on each tumor sample run individually (against the single normal sample), this script will load the resulting data for each sample and merge them, and run multi-sample segmentation via ASCAT's multipcf algorithm. It will then format the output for use with CNAlignR.
##'
##' @export
get_CNAlignR_obj_for_snp_data <- function(ascat_dir, sex, build, normal_sample, GCcontentfile=NULL, replictimingfile=NULL, multipcf_penalty=300, multipcf_refine=F, multipcf_selectAlg='exact', cleanup=T, seed=NA, output_dir='.', obj_filename='CNAlignR_obj.Rdata', Tumor_LogR_filename='Tumor_LogR.txt', Tumor_BAF_filename='Tumor_BAF.txt', Germline_LogR_filename='Germline_LogR.txt', Germline_BAF_filename='Germline_BAF.txt') {

    ascat_files <- dir(ascat_dir,full.names=T)
    Tumor_LogR_files <- grep(paste0(normal_sample,'_Tumor_LogR.txt'), ascat_files, value=T)
    Tumor_BAF_files <- grep(paste0(normal_sample,'_Tumor_BAF.txt'), ascat_files, value=T)
    Germline_LogR_files <- grep(paste0(normal_sample,'_Germline_LogR.txt'), ascat_files, value=T)
    Germline_BAF_files <- grep(paste0(normal_sample,'_Germline_BAF.txt'), ascat_files, value=T)

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

    ## get LogR data 
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

    ## get BAF data 
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

    ## subset for positions in both logR and BAF data
    all_samples <- names(x_LogR)[!names(x_LogR) %in% c('V1','Chromosome','Position')]
    tumor_samples <- all_samples[!all_samples %in% normal_sample]
    valid_rows <- intersect(rownames(x_BAF), rownames(x_LogR))
    x_LogR <- x_LogR[valid_rows,]
    x_BAF <- x_BAF[valid_rows,]
    t_LogR <- x_LogR[,c('Chromosome','Position',tumor_samples)]
    t_BAF <- x_BAF[,c('Chromosome','Position',tumor_samples)]
    n_LogR <- x_LogR[,c('Chromosome','Position',normal_sample)]
    n_BAF <- x_BAF[,c('Chromosome','Position',normal_sample)]

    if(!dir.exists(output_dir)) dir.create(output_dir)
    Merged_Tumor_LogR_filename <- file.path(output_dir,Tumor_LogR_filename)
    Merged_Tumor_BAF_filename <- file.path(output_dir,Tumor_BAF_filename)
    Merged_Germline_LogR_filename <- file.path(output_dir,Germline_LogR_filename)
    Merged_Germline_BAF_filename <- file.path(output_dir,Germline_BAF_filename)
    write.table(t_LogR, file = Merged_Tumor_LogR_filename, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_LogR, file = Merged_Germline_LogR_filename, sep = "\t", quote = FALSE, col.names = NA)
    write.table(t_BAF, file = Merged_Tumor_BAF_filename, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_BAF, file = Merged_Germline_BAF_filename, sep = "\t", quote = FALSE, col.names = NA)

    ## get snp-level data
    message('Running prep_data_for_multipcf() ...')
    obj <- prep_data_for_multipcf(Tumor_LogR_file=Merged_Tumor_LogR_filename,
                                  Tumor_BAF_file=Merged_Tumor_BAF_filename,
                                  Germline_LogR_file=Merged_Germline_LogR_filename,
                                  Germline_BAF_file=Merged_Germline_BAF_filename,
                                  sex=sex,
                                  build=build,
                                  GCcontentfile=GCcontentfile,
                                  replictimingfile=replictimingfile)

    ## add a list with the parameters to the output object so that we have a record of what values were used
    main_params <- list(ascat_dir=ascat_dir, sex=sex, build=build, normal_sample=normal_sample, GCcontentfile=GCcontentfile, replictimingfile=replictimingfile, output_dir=output_dir, obj_filename=obj_filename, Tumor_LogR_filename=Tumor_LogR_filename, Tumor_BAF_filename=Tumor_BAF_filename, Germline_LogR_filename=Germline_LogR_filename, Germline_BAF_filename=Germline_BAF_filename)

    obj$main_params <- main_params
   
    # save object 
    output_path <- file.path(output_dir, obj_filename)
    message('Saving CNAlignR data object to: ', output_path)
    saveRDS(obj, file=output_path)
    
}
 

