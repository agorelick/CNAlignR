##' save_preprocessed_bin_data_for_ascat
##'
##' Write preprocessed bin data into the format expected for ascat.loadData()
##'
##' @export
save_preprocessed_bin_data_for_ascat <- function(preprocessed_data, normal_sample, sex, Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file, Germline_BAF_file) {

    if(sex=='XX') {
        valid_chr <- c(1:22,'X')
    } else {
        valid_chr <- c(1:22,'X','Y')
    }
    d <- preprocessed_data[Chromosome %in% valid_chr,]
    d$Chromosome <- factor(d$Chromosome, levels=valid_chr)
    tumor_samples <- unique(d[sample!=normal_sample,(sample)])

    message('Saving preprocessed bin data for input to ASCAT.')
    t_LogR <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='LogR')
    t_LogR <- t_LogR[order(Chromosome, Position),]
    n_LogR <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='LogR')
    n_LogR <- n_LogR[order(Chromosome, Position),]
    t_BAF <- data.table::dcast(Chromosome + Position ~ sample, data=d[sample %in% tumor_samples], value.var='BAF')
    t_BAF <- t_BAF[order(Chromosome, Position),]
    n_BAF <- dcast(Chromosome + Position ~ sample, data=d[sample %in% normal_sample,], value.var='BAF')
    n_BAF <- n_BAF[order(Chromosome, Position),]

    write.table(t_LogR, file = Tumor_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_LogR, file = Germline_LogR_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(t_BAF, file = Tumor_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
    write.table(n_BAF, file = Germline_BAF_file, sep = "\t", quote = FALSE, col.names = NA)
}
