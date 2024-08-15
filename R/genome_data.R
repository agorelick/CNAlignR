##' genome_data
##' load internal data for b37 or hg38 and return formatted data useful making genome plots and analyses
##' @export
genome_data <- function(build) { 
    if(build=='hg19') {
        data(chr_lengths_b37, envir=environment())
        chr <- chr_lengths_b37
        rm(chr_lengths_b37)
        data(cytoBand_b37, envir=environment())
        bands <- cytoBand_b37
        rm(cytoBand_b37) 
    } else if(build=='hg38') {
        data(chr_lengths_hg38, envir=environment())
        chr <- chr_lengths_hg38
        rm(chr_lengths_hg38)
        data(cytoBand_hg38, envir=environment())
        bands <- cytoBand_hg38
        rm(cytoBand_hg38) 
    }

    ## get a table of hg19 chromosome lengths with global start/end/midpoint positions (in Mb) 
    chr$length <- chr$length / 1e6
    chr$chr_start <- 0
    chr$chr_end <- chr$length
    chr$global_start <- rep(as.integer(NA),nrow(chr))
    chr$global_end <- rep(as.integer(NA),nrow(chr))
    chr$global_midpoint <- rep(as.integer(NA),nrow(chr))
    chr$global_start[1] <- 0
    chr$global_end[1] <- chr$length[1]
    chr$global_midpoint[1] <- (chr$chr_start[1] + chr$chr_end[1])/2
    for(i in 2:nrow(chr)) {
        chr$global_start[i] <- sum(chr$length[1:(i-1)])
        chr$global_end[i] <- chr$global_start[i]+chr$length[i]
        chr$global_midpoint[i] <- (chr$global_start[i] + chr$global_end[i])/2
    }
    chr$chr <- factor(chr$chr, levels=c(1:22,'X','Y','MT'))
    setkey(chr,'chr','chr_start','chr_end')

    names(bands) <- c('chr','band_start','band_end','region','stain')
    bands[,arm:=strtrim(region,1)]
    bands[,chr:=gsub('chr','',chr)]
    bands[,c('stain','region'):=NULL]
    collapse_arm <- function(bands) {
        arm_start <- min(bands$band_start) / 1e6
        arm_end <- max(bands$band_end) / 1e6
        list(arm_start=arm_start, arm_end=arm_end) 
    }
    arms <- bands[,collapse_arm(.SD),by=c('chr','arm')]
    arms$chr <- factor(arms$chr, levels=c(1:22,'X','Y','MT'))    
    arms <- merge(arms, chr, by='chr', all=T)
    arms[arm=='p', arm_start:=chr_start] 
    arms[arm=='q', arm_end:=chr_end] 
    #arms[,global_arm_start:=global_start + arm_start]
    #arms[,global_arm_end:=global_start + arm_end]
    arms <- arms[,c('chr','arm','arm_start','arm_end'),with=F]
    arms[chr=='MT', arm:='']
    arms[chr=='MT',arm_start:=0.0]
    arms[chr=='MT',arm_end:=16569/1e6]
    arms <- arms[!is.na(chr),]
    setkey(arms,'chr','arm_start','arm_end')
    list(chr=chr, arms=arms)
}



