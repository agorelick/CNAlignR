prep_data_for_multipcf <-
function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file, Germline_BAF_file, sex, genomeVersion, GCcontentfile=NULL, replictimingfile=NULL, build='hg19') {

    require(ASCAT)
    n_LogR <- fread(Germline_LogR_file)
    if(sex=='XX') {
        gender = rep('XX',nrow(n_LogR))
    } else if(sex=='XY') {
        gender = rep('XY',nrow(n_LogR))
    } else {
        gender <- NULL
    }
    chr_levels <- c(1:22,'X','Y')

    message('Loading ASCAT data ...')
    ascat.bc = ascat.loadData(Tumor_LogR_file = Tumor_LogR_file,
                              Tumor_BAF_file = Tumor_BAF_file,
                              Germline_LogR_file = Germline_LogR_file,
                              Germline_BAF_file = Germline_BAF_file,
                              gender = gender, 
                              genomeVersion = build) 
    ascat.loadData.params <- list(Tumor_LogR_file=Tumor_LogR_file, Tumor_BAF_file=Tumor_BAF_file, Germline_LogR_file=Germline_LogR_file, Germline_BAF_file=Germline_BAF_file, sex=sex, genomeVersion=genomeVersion, gender = unique(gender), genomeVersion = build)
   
    ## get list of tumor samples and (single) normal sample used throughout the pipeline
    tumor_samples <- names(ascat.bc$Tumor_LogR)
    normal_sample <- names(ascat.bc$Germline_LogR)

    ## correct LogR for GC content and replication timing
    if(!is.null(GCcontentfile) & !is.null(replictimingfile)) {
        message('Correcting LogR (may take a few minutes) ...')
        ascat.bc.corrected = ascat.correctLogR(ascat.bc, GCcontentfile = GCcontentfile,  replictimingfile=replictimingfile)
    } else {
        ascat.bc.corrected <- ascat.bc
    }
    ascat.correctLogR.params <- list(GCcontentfile=GCcontentfile, replictimingfile=replictimingfile)

    ## format data from ASCAT
    message('Processing marker-level data ...')

    ## get chr-arm positions and SNPpos data for annotating markers
    arms <- genome_data(build)$arms 
    arms[,arm_start:=arm_start*1e6]
    arms[,arm_end:=arm_end*1e6]
    setkey(arms,'chr','arm_start','arm_end')
    SNPpos <- ascat.bc.corrected$SNPpos

    ## annotate the LogR data
    t_LogR <- ascat.bc.corrected$Tumor_LogR
    n_LogR <- ascat.bc.corrected$Germline_LogR
    dt_LogR <- cbind(SNPpos, t_LogR, n_LogR) 
    dt_LogR <- cbind(VariantID=rownames(dt_LogR), as.data.table(dt_LogR))
    dt_LogR[,Position2:=Position+1]
    dt_LogR$Chromosome <- factor(dt_LogR$Chromosome, levels=chr_levels)
    setkey(dt_LogR,'Chromosome','Position','Position2')
    dt_LogR <- foverlaps(dt_LogR, arms, type='within')
    dt_LogR <- dt_LogR[order(Chromosome, Position),]
    dt_LogR[,c('Position2','arm_start','arm_end'):=NULL]
    dt_LogR <- data.table::melt(dt_LogR, id.vars=c('Chromosome','arm','Position','VariantID'))
    setnames(dt_LogR,c('variable','value'),c('sample','LogR'))

    ## annotate the BAF data
    t_BAF <- ascat.bc.corrected$Tumor_BAF
    n_BAF <- ascat.bc.corrected$Germline_BAF
    dt_BAF <- cbind(SNPpos, t_BAF, n_BAF) 
    dt_BAF <- cbind(VariantID=rownames(dt_BAF), as.data.table(dt_BAF))
    dt_BAF[,Position2:=Position+1]
    dt_BAF$Chromosome <- factor(dt_BAF$Chromosome, levels=chr_levels)
    setkey(dt_BAF,'Chromosome','Position','Position2')
    dt_BAF <- foverlaps(dt_BAF, arms, type='within')
    dt_BAF <- dt_BAF[order(Chromosome, Position),]
    dt_BAF[,c('Position2','arm_start','arm_end'):=NULL]
    dt_BAF <- data.table::melt(dt_BAF, id.vars=c('Chromosome','arm','Position','VariantID'))
    setnames(dt_BAF,c('variable','value'),c('sample','BAF'))

    ## combine LogR and BAF data
    dt <- cbind(dt_LogR, BAF=dt_BAF$BAF)  # need to add error-check to ensure data lines up (unless I use inefficient merging below)
    #dt <- merge(dt_LogR, dt_BAF[,c('sample','VariantID','BAF'),with=F], by=c('sample','VariantID'), all=T)
    dt[Chromosome=='MT', arm:='']
    dt <- dt[!is.na(arm),]
    dt[,charm:=paste0(Chromosome,arm)]
    dt <- dt[order(Chromosome, Position, sample),]
    dt[is.nan(BAF),BAF:=NA]
    dt[is.na(BAF),Ref:='.']
    dt[is.na(BAF),Alt:='.']
    dt[!is.na(BAF),Ref:='N']
    dt[!is.na(BAF),Alt:='N']
    dt[,count1:=as.numeric(NA)]
    dt[,count2:=as.numeric(NA)]
    dt[,count_raw:=as.numeric(NA)]
    dt[,dp:=as.numeric(NA)]
    dt <- dt[,c('Chromosome','arm','Position','Ref','Alt','sample','count1','count2','count_raw','LogR','dp','BAF','VariantID','charm'),with=F]
    dt$Chromosome <- factor(dt$Chromosome, levels=chr_levels)
    dt <- dt[order(Chromosome,Position,sample),]
    dt <- dt[!sample %in% normal_sample,]

    ## get number of germline copies (germline_copies) for each chromosome including sex chromosomes
    if(sex=='XX') {
        dt[Chromosome %in% 'X',germline_copies:=2]
        dt[Chromosome %in% 'Y',germline_copies:=0]
    } else if(sex=='XY') {
        dt[Chromosome %in% 'X',germline_copies:=1]
        dt[Chromosome %in% 'Y',germline_copies:=1]   
    } else {
        dt[Chromosome %in% c('X','Y'),germline_copies:=NA]
    }
    
    l <- split(dt, by='sample')

    message('Done!')
    list(marker_level=l, ascat.bc=ascat.bc.corrected, ascat.loadData.params=ascat.loadData.params, ascat.correctLogR.params=ascat.correctLogR.params, tumor_samples=tumor_samples, normal_sample=normal_sample)
}
