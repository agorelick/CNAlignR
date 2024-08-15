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



