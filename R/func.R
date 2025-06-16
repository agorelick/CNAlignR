maf_fields <- function() {
    protected_fields <- c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'One_Consequence', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'TRANSCRIPT_STRAND', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'ExAC_AF', 'ExAC_AF_Adj', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'CONTEXT', 'src_vcf_id', 'tumor_bam_uuid', 'normal_bam_uuid', 'case_id', 'GDC_FILTER', 'COSMIC', 'MC3_Overlap', 'GDC_Validation_Status', 'GDC_Valid_Somatic', 'vcf_region', 'vcf_info', 'vcf_format', 'vcf_tumor_gt', 'vcf_normal_gt')
    protected_fields
}


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


##' get_fits_for_solution
##' @export
get_fits_for_solution <- function(model, sol_number) {
    sol <- paste0('Solution_',sol_number)
    df <- model$df[,c('Variable',sol)]
    df$info <- sub(".*\\[(.*?)\\].*", "\\1", as.character(df$Variable)) 
    pu <- df[grepl('pu\\[',df$Variable),]
    pu$Variable <- NULL
    setnames(pu,sol,'pu')
    pl <- df[grepl('pl\\[',df$Variable),]
    pl$Variable <- NULL
    setnames(pl,sol,'pl')
    fits <- merge(pu, pl, by='info')
    setnames(fits,'info','sample')
    fits
}


##' get_highlighted_segments_for_solution
##' @export
get_highlighted_segments_for_solution <- function(model, sol_number) {
    sol <- paste0('Solution_',sol_number)
    df <- model$df[,c('Variable',sol)]
    df$info <- sub(".*\\[(.*?)\\].*", "\\1", as.character(df$Variable)) 
    allmatch <- df[grepl('allmatch\\[',df$Variable),]
    allmatch$Variable <- NULL
    setnames(allmatch,sol,'highlight')
    highlighted <- allmatch$info[allmatch$highlight==1]
    highlighted <- as.integer(gsub('seg','',highlighted))
    highlighted
}


##' get_fit_and_plot
##'
##' Plot segments for a given sample after CNAlign
##'
##' @export
get_fit_and_plot <- function(s, highlight, fits, obj, params, plot_dir) { 
    # create a directory for the parameters
    params <- params[!names(params) %in% c('gurobi_license','normal_baseline')]
    params$rho <- round(params$rho, 3)
    #browser()

    names(params)[names(params)=='rho'] = 'r'
    names(params)[names(params)=='delta_mcn_to_avg'] = 'dm2ma'
    names(params)[names(params)=='delta_mcn_to_int'] = 'dm2i'
    names(params)[names(params)=='delta_mcnavg_to_int'] = 'dma2i'
    names(params)[names(params)=='delta_tcn_to_avg'] = 'dt2ta'
    names(params)[names(params)=='delta_tcn_to_int'] = 'dt2i'
    names(params)[names(params)=='delta_tcnavg_to_int'] = 'dta2i'
    names(params)[names(params)=='max_homdel_mb'] = 'homdel'
    names(params)[names(params)=='min_ploidy'] = 'pl1'
    names(params)[names(params)=='max_ploidy'] = 'pl2'
    names(params)[names(params)=='min_purity'] = 'pu1'
    names(params)[names(params)=='max_purity'] = 'pu2'
    names(params)[names(params)=='mcn_weight'] = 'mw'
    names(params)[names(params)=='min_aligned_seg_mb'] = 'minlen'
    names(params)[names(params)=='min_cna_segments_per_sample'] = 'mincna'
    names(params)[names(params)=='timeout'] = 'sec'
    names(params)[names(params)=='obj2_clonalonlyonly'] = 'clonal'

    outdir <- file.path(plot_dir,paste(paste0(names(params),'=',params), collapse='_'))
    if(!dir.exists(outdir)) dir.create(outdir)

    # figure suffix has purity/ploidy
    fits <- as.data.table(fits)
    pu <- round(fits[sample==s,(pu)], 3)
    pl <- round(fits[sample==s,(pl)], 3)
    suffix <- paste0('pu=',pu,'_pl=',pl,'.png')

    message(s,': pu=',pu,', pl=',pl)
    fit <- get_fit(s, obj, purity=pu, ploidy=pl)
    p <- plot_fit(s, fit, obj, highlight_seg=highlight, LogR_point_size=0.25, BAF_point_size=0.25)

    ggsave(plot=p, file.path(outdir,paste(s,suffix,sep='_')),width=8,height=8)
    fit
}


##' plot_tcn_angular_distance_heatmap
##'
##' Plot a preliminary copy number heatmap based on angular distances from total copy number
##'
##' @export
plot_tcn_angular_distance_heatmap <- function(obj=NULL, segs=NULL, build=NULL, groups=NULL, group_colors=NULL, normal_sample, patient='') {

    if(!is.null(obj)) {
        segs <- rbindlist(obj$segment_level, fill=T)
        build <- obj$main_params$build
        patient <- obj$main_params$patient
        patient <- obj$main_params$patient
    } 
    
    sex <- obj$main_params$sex
    if(sex=='XX') {
        valid_chr <- c(1:22,'X')
    } else {
        valid_chr <- c(1:22,'X','Y')
    }

    if(is.null(obj) & (is.null(segs) | is.null(build))) stop('Must provide either [obj] OR ([segs] AND [build])') 

    ## angular distance tree from segment LogRs
    LogR_mat <- d2m(data.table::dcast(segment ~ sample, value.var='LogR_segmented', data=segs))
    LogR_mat <- LogR_mat[!is.na(rowSums(LogR_mat)),]
    n_seg <- nrow(LogR_mat)
    n_samp <- ncol(LogR_mat)
    for(i in 1:n_samp) LogR_mat[,i] <- LogR_mat[,i] / sqrt(sum(LogR_mat[,i]^2))

    long <- as.data.table(reshape2::melt(LogR_mat))
    names(long) <- c('segment','sample','value')
    segs <- merge(long, obj$segments, by='segment', all.x=T)
    segs <- segs[Chromosome %in% valid_chr,]
    segs$Chromosome <- factor(segs$Chromosome, levels=valid_chr)

    ## get angular distance matrix/tree
    numerator <- t(LogR_mat) %*% LogR_mat
    magnitudes <- sqrt(colSums(LogR_mat^2))
    denominator <- magnitudes %*% t(magnitudes)
    cos_theta_LogR <- numerator / denominator
    theta_LogR <- acos(cos_theta_LogR)
    for(i in 1:ncol(theta_LogR)) theta_LogR[i,i] <- 0
    theta_LogR <- cbind(theta_LogR, diploid=pi/3)
    theta_LogR <- rbind(theta_LogR, diploid=pi/3)
    for(i in 1:(ncol(LogR_mat)+1)) theta_LogR[i,i] <- 0
    tree_LogR <- nj(theta_LogR)
    tree_LogR <- phytools::reroot(tree_LogR, which(tree_LogR$tip.label=='diploid'))
    tree_LogR$tip.label[tree_LogR$tip.label=='diploid'] <- normal_sample

    setnames(segs,c('Chromosome'),c('chr'))
    segs$sample <- as.character(segs$sample)
    message('Expanding segments to include NA regions in each chromosome ...')

    expand_segments_to_complete_chromosome_for_sample <- function(this.sample, segs, gr_chr) {
        #browser()
        message(this.sample)
        mat_sample <- segs[sample==this.sample,]
        mat_sample$chr <- factor(mat_sample$chr, levels(seqnames(gr_chr)))
        mat_sample <- mat_sample[order(chr)]
        gr_mat_sample <- makeGRangesFromDataFrame(mat_sample,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='seg_start', end.field='seg_end')
        NA_regions <- BiocGenerics::setdiff(gr_chr, gr_mat_sample)
        complete_regions <- as.data.table(sort(c(gr_mat_sample, NA_regions)))
        complete_regions$sample <- this.sample
        complete_regions[,segment:=paste0(seqnames,':',start,'-',end)]
        setnames(complete_regions,'seqnames','chr')
        complete_regions <- complete_regions[order(chr, start, end),]
        ## collapse regions with no difference in copy number
        complete_regions$sample <- this.sample
        complete_regions
    }

    gd <- genome_data(build)$chr
    gd <- gd[chr %in% valid_chr]
    gd$chr <- factor(gd$chr, levels=valid_chr)
    gd[,chr_start:=chr_start * 1e6]
    gd[,chr_end:=chr_end * 1e6]

    gr_chr <- makeGRangesFromDataFrame(gd,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='chr_start', end.field='chr_end')
    sample_list <- lapply(unique(segs$sample), expand_segments_to_complete_chromosome_for_sample, segs, gr_chr)
    segs <- rbindlist(sample_list)

    segs <- merge(segs, gd[,c('chr','global_start'),with=F], by='chr', all.x=T)
    segs[,global_seg_start_mb:=global_start + start/1e6]
    segs[,global_seg_end_mb:=global_start + end/1e6]
    segs <- segs[chr %in% valid_chr]
    segs$chr <- factor(segs$chr, levels=valid_chr)
    segs[,segment:=paste0(chr,':',start,'-',end)]

    ## add values for the diploid normal
    toadd <- segs[!duplicated(segment),]
    toadd[,value:=0]
    toadd[,sample:=normal_sample]
    segs <- rbind(segs, toadd)

    p1 <- ggtree(tree_LogR, linewidth=0.5)
    p1 <- p1 + guides(fill='none') + theme(legend.position='none') + labs(title=patient, subtitle='CN angular distance tree')
    pd <- as.data.table(p1$data)
    pd <- pd[isTip==T,]
    pd <- pd[order(y,decreasing=F),]
    pd$label <- factor(pd$label, levels=pd$label)
    segs$sample <- factor(segs$sample, levels=pd$label)
    segs$sample.i <- as.integer(segs$sample)

    if(!is.null(groups)) {
        message('Adding sample type groups')
        p1 <- p1 %<+% groups[,c('lesion_id','group'),with=F]
        p1 <- p1 + geom_tiplab(aes(color=group))
        p1 <- p1 + scale_color_manual(values=group_colors,name='Organ category')
    } else {
        p1 <- p1 + geom_tiplab()
    }

    p1 <- p1 + xlim(0, 1.15*max(p1$data$x))
    gd2 <- copy(gd)
    chr19_22 <- gd[chr==19,]
    chr19_22$length <- sum(gd[chr %in% c(19,20,21,22),(length)])
    chr19_22[,global_midpoint:=global_start+0.5*length]
    chr19_22[,global_end:=global_start+length]
    chr19_22$chr <- '19-22'
    gd2 <- gd2[!chr %in% c(19,20,21,22)]
    gd2 <- rbind(gd2, chr19_22) 
    gd2 <- gd2[order(global_midpoint),]

    right_labs <- data.table(sample=levels(segs$sample))
    right_labs$sample <- factor(right_labs$sample, levels=right_labs$sample)
    right_labs$sample.i <- as.integer(right_labs$sample)
    right_labs[,y:=sample.i]

    p2 <- ggplot(segs) + 
        scale_x_continuous(expand=c(0,-0.1), breaks=gd2$global_midpoint, labels=gd2$chr) + 
        scale_y_continuous(expand=c(0,0), breaks=right_labs$y, labels=right_labs$sample, position='right') + 
        geom_rect(aes(ymin=sample.i-0.5, ymax=sample.i+0.5, xmin=global_seg_start_mb, xmax=global_seg_end_mb, fill=value)) + 
        scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, na.value='#bfbfbf',name='Normalized TCN') +
        geom_vline(xintercept=c(0,tail(gd$global_end,1)), linewidth=0.5, linetype='solid') +
        geom_hline(yintercept=seq(min(segs$sample.i)-0.5,max(segs$sample.i)+0.5), linewidth=0.25, linetype='solid') + 
        geom_vline(xintercept=tail(gd$global_start,-1), linewidth=0.25, linetype='dotted') + 
        guides(color='none') +
        theme_ang(base_size=10) +
        labs(x='Genomic Position', y=NULL, title='', subtitle='Normalized total copy number') +
        theme(axis.line=element_blank(), legend.position='bottom', axis.ticks=element_blank(), axis.text.y=element_blank())

    p <- plot_grid(p1, p2, nrow=1, rel_widths=c(1.5,3), align='h', axis='tb')
    p
}



##' get_ascn_segments
##'
##' Extract a table of copy number segments with purity/ploidy-corrected allele-specific copy number
##'
##' @export
get_ascn_segments <- function(obj, fit_file, normal_sample, min_purity=0.1) {
    message('Extracting purity/ploidy-corrected ASCN segment data ...')
    fits <- fread(fit_file)
    fits <- fits[pu >= min_purity,]
    segs <- rbindlist(obj$segment_level, fill=T)
    get_fits <- function(i, fits, segs) {
        pu <- fits$pu[i]
        pl <- fits$pl[i]
        this.sample <- fits$sample[i]
        message(this.sample)
        x <- segs[sample==this.sample]
        out <- get_na_nb(pu, pl, x, rfield='LogR_segmented', bfield='BAF_segmented')
        out$pu <- pu
        out$pl <- pl
        out
    }
    fits <- fits[order(sample),]
    l <- lapply(1:nrow(fits), get_fits, fits, segs)
    segs <- rbindlist(l)

    ## add Normal segment pseudo-values
    toadd <- segs[!duplicated(segment),]
    toadd[,sample:=normal_sample]
    toadd[Chromosome %in% 1:22,c('na','nb'):=list(1,1)]
    sex <- obj$main_params$sex
    if(sex=='XX') {
        toadd[Chromosome=='X',c('na','nb'):=list(1,1)]
        toadd[Chromosome=='Y',c('na','nb'):=list(0,0)]
    } else if(sex=='XY') {
        toadd[Chromosome=='X',c('na','nb'):=list(1,0)]
        toadd[Chromosome=='Y',c('na','nb'):=list(1,0)]
    } else if(is.null(sex)) {
        toadd[Chromosome=='X',c('na','nb'):=list(NA,NA)]
        toadd[Chromosome=='Y',c('na','nb'):=list(NA,NA)]   
    }
    toadd[,LogR_segmented:=NA]
    toadd[,BAF_segmented:=NA]
    toadd$pu <- 0; toadd$pl <- 2
    segs <- rbind(segs, toadd)
    segs[!is.na(na) & !is.na(nb),tcn:=na+nb]
    segs[!is.na(na) & is.na(nb),tcn:=na]
    segs
}


##' plot_tcn_heatmap
##'
##' generate a heatmap with total copy number, a NJ tree based on TCN, and a barplot annotating tumor purity
##'
##' @export
plot_tcn_heatmap <- function(segs, fits, groups, group_colors, sex, build, normal_sample) {

    ## angular distance tree from segment LogRs
    tcn_mat <- d2m(data.table::dcast(segment ~ sample, value.var='tcn', data=segs))
    tcn_mat <- tcn_mat[!is.na(rowSums(tcn_mat)),]

    long <- as.data.table(reshape2::melt(tcn_mat))
    names(long) <- c('segment','sample','value')
    segs <- merge(long, obj$segments, by='segment', all.x=T)

    ## get angular distance matrix/tree
    dm <- dist(t(tcn_mat), method='manhattan')
    tree_tcn <- nj(dm)
    tree_tcn <- phytools::reroot(tree_tcn, which(tree_tcn$tip.label==normal_sample))

    setnames(segs,c('Chromosome'),c('chr'))
    segs$sample <- as.character(segs$sample)
    message('Expanding segments to include NA regions in each chromosome ...')

    if(sex=='XX') {
        valid_chr <- c(1:22,'X')
    } else {
        valid_chr <- c(1:22,'X','Y')
    }

    expand_segments_to_complete_chromosome_for_sample <- function(this.sample, segs, gr_chr) {
        message(this.sample)
        mat_sample <- segs[sample==this.sample,]
        mat_sample$chr <- factor(mat_sample$chr, levels(seqnames(gr_chr)))
        mat_sample <- mat_sample[order(chr)]
        gr_mat_sample <- makeGRangesFromDataFrame(mat_sample,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='seg_start', end.field='seg_end')
        NA_regions <- BiocGenerics::setdiff(gr_chr, gr_mat_sample)
        complete_regions <- as.data.table(sort(c(gr_mat_sample, NA_regions)))
        complete_regions$sample <- this.sample
        complete_regions[,segment:=paste0(seqnames,':',start,'-',end)]
        setnames(complete_regions,'seqnames','chr')
        complete_regions <- complete_regions[order(chr, start, end),]
        ## collapse regions with no difference in copy number
        complete_regions$sample <- this.sample
        complete_regions
    }

    gd <- genome_data(build)$chr
    gd <- gd[chr %in% valid_chr]
    gd$chr <- factor(gd$chr, levels=valid_chr)
    gd[,chr_start:=chr_start * 1e6]
    gd[,chr_end:=chr_end * 1e6]

    segs <- segs[chr %in% valid_chr,]
    segs$chr <- factor(segs$chr, levels=valid_chr)
    gr_chr <- makeGRangesFromDataFrame(gd,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='chr_start', end.field='chr_end')
    sample_list <- lapply(unique(segs$sample), expand_segments_to_complete_chromosome_for_sample, segs, gr_chr)
    segs <- rbindlist(sample_list)

    segs <- merge(segs, gd[,c('chr','global_start'),with=F], by='chr', all.x=T)
    segs[,global_seg_start_mb:=global_start + start/1e6]
    segs[,global_seg_end_mb:=global_start + end/1e6]
    segs <- segs[chr %in% valid_chr]
    segs$chr <- factor(segs$chr, levels=valid_chr)
    segs[,segment:=paste0(chr,':',start,'-',end)]

    blues <- brewer.pal(9,'RdBu')[c(9,7)]
    reds <-  c(brewer.pal(9, 'OrRd')[c(1,3,5,7,9)],'black')
    cols <- c(blues, 'white', reds,'#bfbfbf')
    names(cols) <- c(seq(0,7,by=1),'8+','n/a')

    segs$tcn.i <- round(segs$value)
    segs$tcn.f <- as.character(segs$tcn.i)
    segs[is.na(tcn.i), tcn.f:='n/a']
    segs[tcn.i >= 8, tcn.f:='8+']
    segs[tcn.i < 0, tcn.f:='0']
    segs$tcn.f <- factor(segs$tcn.f, levels=c(0:7,'8+','n/a'))
    segs[value > 8, value:=8]

    p1 <- ggtree(tree_tcn, linewidth=0.5)
    p1 <- p1 + guides(fill='none') + theme(legend.position='none') + labs(title=patient, subtitle='TCN NJ tree')
    p1 <- p1 %<+% groups
    pd <- as.data.table(p1$data)
    pd <- pd[isTip==T,]
    pd <- pd[order(y,decreasing=F),]
    pd$label <- factor(pd$label, levels=pd$label)
    segs$sample <- factor(segs$sample, levels=pd$label)
    segs$sample.i <- as.integer(segs$sample)
    p1 <- p1 + geom_tiplab(aes(color=group))
    p1 <- p1 + scale_color_manual(values=group_colors,name='Sample type')
    p1 <- p1 + xlim(0, 1.15*max(p1$data$x))

    gd2 <- copy(gd)
    chr19_22 <- gd[chr==19,]
    chr19_22$length <- sum(gd[chr %in% c(19,20,21,22),(length)])
    chr19_22[,global_midpoint:=global_start+0.5*length]
    chr19_22[,global_end:=global_start+length]
    chr19_22$chr <- '19-22'
    gd2 <- gd2[!chr %in% c(19,20,21,22)]
    gd2 <- rbind(gd2, chr19_22) 
    gd2 <- gd2[order(global_midpoint),]

    right_labs <- data.table(sample=levels(segs$sample))
    right_labs$sample <- factor(right_labs$sample, levels=right_labs$sample)
    right_labs$sample.i <- as.integer(right_labs$sample)
    right_labs[,y:=sample.i]

    p2 <- ggplot(segs) + 
        scale_x_continuous(expand=c(0,-0.1), breaks=gd2$global_midpoint, labels=gd2$chr) + 
        scale_y_continuous(expand=c(0,0), breaks=right_labs$y, labels=right_labs$sample, position='right') + 
        geom_rect(aes(ymin=sample.i-0.5, ymax=sample.i+0.5, xmin=global_seg_start_mb, xmax=global_seg_end_mb, fill=tcn.f)) + 
        scale_fill_manual(values=cols, na.value='#bfbfbf',name='TCN (nearest integer)') +
        #scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=2, na.value='#bfbfbf',name='True TCN') +
        geom_vline(xintercept=c(0,tail(gd$global_end,1)), linewidth=0.5, linetype='solid') +
        geom_hline(yintercept=seq(min(segs$sample.i)-0.5,max(segs$sample.i)+0.5), linewidth=0.25, linetype='solid') + 
        geom_vline(xintercept=tail(gd$global_start,-1), linewidth=0.25, linetype='dotted') + 
        guides(color='none') +
        theme_ang(base_size=10) +
        labs(x='Genomic Position', y=NULL, title='', subtitle='Purity/ploidy-corrected total copy number') +
        theme(axis.line=element_blank(), legend.position='bottom', axis.ticks=element_blank(), axis.text.y=element_blank())

    pd_fits <- data.table(sample=right_labs$sample)
    pd_fits <- merge(pd_fits, fits, by='sample', all.x=T)
    pd_fits[sample==normal_sample,c('pl','pu'):=list(2,0)]
    pd_fits$sample <- factor(pd_fits$sample, levels=right_labs$sample)
    p3 <- ggplot(pd_fits, aes(x=sample, y=pu)) +  
        scale_y_continuous(expand=c(0,0)) + 
        geom_bar(stat='identity')  +
        coord_flip() +
        theme_ang(base_size=10) +
        labs(x=NULL, y='Purity') + 
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

    p <- plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1.5,3,0.5), align='h', axis='tb')
    p
}




