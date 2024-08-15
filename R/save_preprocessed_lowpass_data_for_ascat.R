##' save_preprocessed_lowpass_data_for_ascat
##' Write preprocessed lowpass data into the format expected for ascat.loadData()
##' @export
save_preprocessed_lowpass_data_for_ascat <- function(d, normal_sample, sex, tmpdir='.') { 
    require(ASCAT)

    if(sex=='XX') {
        valid_chr <- c(1:22,'X')
    } else {
        valid_chr <- c(1:22,'X','Y')
    }
    d <- d[Chromosome %in% valid_chr,]
    d$Chromosome <- factor(d$Chromosome, levels=valid_chr)
    tumor_samples <- unique(d[sample!=normal_sample,(sample)])

    message('Saving preprocessed data from lpASCN for input to ASCAT.')
    t_LogR <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='LogR')
    t_LogR <- t_LogR[order(Chromosome, Position),]
    n_LogR <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='LogR')
    n_LogR <- n_LogR[order(Chromosome, Position),]
    t_BAF <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='BAF')
    t_BAF <- t_BAF[order(Chromosome, Position),]
    n_BAF <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='BAF')
    n_BAF <- n_BAF[order(Chromosome, Position),]

    if(!dir.exists(tmpdir)) dir.create(tmpdir)
    tumor_LogR_file <- file.path(tmpdir,'Tumor_LogR.txt')
    tumor_BAF_file <- file.path(tmpdir,'Tumor_BAF.txt')
    germline_LogR_file <- file.path(tmpdir,'Germline_LogR.txt')
    germline_BAF_file <- file.path(tmpdir,'Germline_BAF.txt')
    write.table(t_LogR, file = tumor_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_LogR, file = germline_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(t_BAF, file = tumor_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_BAF, file = germline_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
}
