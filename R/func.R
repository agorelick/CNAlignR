##' theme_ang
##' @export
theme_ang <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
        theme(line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"), 
              text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = F), 
              axis.text = element_text(colour = "black", size = rel(0.8)), 
              axis.ticks = element_line(colour = "black", linewidth = rel(1)), 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = rel(1)),
              legend.key = element_blank(), 
              strip.background = element_blank())
}


##' block_phased_snps
##' @export
block_phased_snps <- function(phased_snps, max_phaseable_distance) {
    ## within each bin, we create phasing blocks 
    pos <- phased_snps$POS
    bin.i <- phased_snps$bin.i
    blocks <- phaseblock(pos, bin.i, max_phaseable_distance)
    phased_snps$block <- blocks
    phased_snps
}


##' dir.create.new
##' @export
dir.create.new <- function(path) {
    if(!dir.exists(path)) {
        message('Creating directory: ', path)
        dir.create(path, recursive=T)
    } else {
        message(path,' already exists.')
    }
}


##' write_distance_matrix
##' @export
write_distance_matrix <- function (dm, file) {
    write.table(dm, file = file, sep = "\t", quote = FALSE, col.names = NA)
}


##' read_distance_matrix
##' @export
read_distance_matrix <- function (file, return.as.matrix = T) {
    distance_matrix <- fread(file)
    rows <- distance_matrix[[1]]
    distance_matrix <- distance_matrix[, (2:ncol(distance_matrix)), with = F]
    m <- as.matrix(distance_matrix)
    rownames(m) <- rows
    if (return.as.matrix == F) {
        as.dist(m, diag = T)
    }
    else {
        m
    }
}

##' write_tsv
##' @export
write_tsv <- function (d, file, sep = "\t", quote = F, row.names = F, ...) {
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}


##' d2m
##' @export
d2m <- function (dt) {
    rows <- dt[[1]]
    if ("data.table" %in% class(dt)) {
        dt <- dt[, c(2:ncol(dt)), with = F]
    }
    else if (class(dt) == "data.frame") {
        dt <- dt[, c(2:ncol(dt))]
    }
    else {
        stop("enter a data.table or a data.frame")
    }
    m <- as.matrix(dt)
    rownames(m) <- rows
    m
}


##' write_refphase_segs
##'
##' Output refphase segments as a tsv-file
##'
##' @export
write_refphase_segs <- function (segs, cn_events = NULL, file, output_format = "summary") {
    columns_to_rename <- c(group_name = "sample_id", seqnames = "chrom")
    d <- as.data.frame(segs)
    d <- S4Vectors::rename(d, columns_to_rename)
    if (output_format == "summary") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "width", "cn_a", "cn_b", "was_cn_updated", "is_ai",
                             "mirrored_vs_ref", "is_reference", "any_ai", "ai_pvalue",
                             "effect_size", "diptest_pvalue", "heterozygous_SNP_number",
                             "homozygous_SNP_number")
        d_filt <- d[, columns_to_save]

    }
    else if (output_format == "full") {
        d <- merge(x = d, y = cn_events, by.x = c("sample_id",
                                                  "chrom", "start", "end", "width", "cn_a", "cn_b",
                                                  "is_ai", "mirrored_vs_ref", "is_LOH"), by.y = c("sample_id",
                                                  "seqnames", "start", "end", "width", "cn_a", "cn_b",
                                                  "is_ai", "mirrored_vs_ref", "is_LOH"))
        columns_to_save <- c("sample_id", "purity", "ploidy",
                             "chrom", "start", "end", "width", "cn_a", "cn_b",
                             "was_cn_updated", "negative_cn_called", "is_ai",
                             "mirrored_vs_ref", "is_reference", "any_ai", "ai_pvalue",
                             "effect_size", "diptest_pvalue", "heterozygous_SNP_number",
                             "homozygous_SNP_number", "amplification_logRthreshold",
                             "gain_logRthreshold", "loss_logRthreshold", "meanlogR",
                             "is_relative_amplification_mean", "is_relative_gain_mean",
                             "is_relative_loss_mean", "is_relative_amplification_ttest",
                             "is_relative_gain_ttest", "is_relative_loss_ttest",
                             "is_absolute_amplification", "is_absolute_gain",
                             "is_absolute_loss", "is_LOH", "is_homozygous_deletion")
        d_filt <- d[, columns_to_save]

    }
    else if (output_format == "copynumbers") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "cn_a", "cn_b")
        d_filt <- d[, columns_to_save]
    }
    else if (output_format == "copynumbers_int") {
        columns_to_save <- c("sample_id", "chrom", "start", "end",
                             "cn_a_integer", "cn_b_integer")
        d_filt <- d[, columns_to_save]
        names(d_filt) <- gsub('_integer','',names(d_filt))
    }
    else {
        stop("output_format has to be in ['summary', 'full', 'copynumbers']")
    }
    utils::write.table(d_filt, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
}



