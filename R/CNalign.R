##' CNalign
##'
##' Determine purity and ploidy values for multiple tumor samples with shared ancestry. This function uses the GuRoBi solver to determine purity/ploidy values for each sample that will *maximize* the number of segments with the same (allele-specific) integer copy numbers in at least rho% of samples. 
##'
##' @export
CNalign <- function(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_homdels=0, max_homdels=3, t_diff=0.05, t_err=0.1, rho=0.85, both_alleles_must_align=1, epsilon=1e-4, assume_wgd=F, tcn_only=F, gurobi_license='~/gurobi.lic', py_script=NA) {
    require(reticulate)
    require(lubridate)
    if(is.na(py_script)) py_script <- system.file("python", "align.py", package = "CNalign")

    ## given input matrices of aligned sample/segment:
    ## 1. logR (required, no NAs)
    ## 2. BAF (optional, can be fully/partially NAs)
    ## Determine purity/ploidy values for each sample such that there is maximum number of segments with the same copy number values for >= r% of the samples

    if(tcn_only==F) {
        source_python(py_script)
        start_time <- now()
        message('Started CNalign at ',as.character(start_time))
        m <- CNalign(dat, min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_homdels=min_homdels, max_homdels=max_homdels, t_diff=t_diff, t_err=t_err, rho=rho, both_alleles_must_align=both_alleles_must_align, gurobi_license=gurobi_license, epsilon=epsilon, assume_wgd=assume_wgd)
    } else {
        source_python(py_script)
        start_time <- now()
        message('Started CNalign at ',as.character(start_time))
        both_alleles_must_align <- NA
        m <- CNalign_tcn(dat, min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_homdels=min_homdels, max_homdels=max_homdels, t_diff=t_diff, t_err=t_err, rho=rho, gurobi_license=gurobi_license, epsilon=epsilon, assume_wgd=assume_wgd)
    }
    end_time <- now()
    message('Ended CNalign at ',as.character(end_time))
    sec_elapsed <- round(as.numeric(end_time) - as.numeric(start_time))
    message('Time elapsed: ',sec_elapsed,'s')

    params <- list(min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_homdels=min_homdels, max_homdels=max_homdels, t_diff=t_diff, t_err=t_err, rho=rho, both_alleles_must_align=both_alleles_must_align, gurobi_license=gurobi_license, epsilon=epsilon, assume_wgd=assume_wgd, tcn_only=tcn_only)

    # return the model object and a list of all the parameters
    out <- list(m=m, params=params)
}


get_solutions <- function(m) {
    get_values_for_solution <- function(solnumber, m) {
        m$setParam('SolutionNumber', solnumber)
        vars <- m$getVars()
        l <- length(vars)
        extract_var <- function(i, vars) {
            name <- vars[[i]]$VarName 
            value <- vars[[i]]$X
            list(name=name, value=value)
        }
        values_for_solution <- rbindlist(lapply(1:l, extract_var, vars))
        values_for_solution$solution <- solnumber
        message('Solution ',solnumber,' extracted.')
        values_for_solution
    }
    solutions <- 0:m$SolCount
    message('Extracting values for ',length(solutions),' solutions ...')
    result_list <- lapply(solutions, get_values_for_solution, m)
    result <- rbindlist(result_list)
    result[,solution:=paste0('sol',solution)]

#    ## floating-point copy number for allele1 per sample per segment
#    parse <- function(name) {
#        n <- strsplit(name,'[,]')[[1]]     
#        list(name=name, sample=n[1], segment=n[2])
#    }
#    n1 <- result[grepl('n1_',name)]
#    n1[,suffix:=gsub('n1_','',name)]
#    n1_parsed <- rbindlist(lapply(unique(n1$name), parse))
#    n1 <- merge(n1, n1_parsed, by='name', all.x=T)
#    n1$segment <- round(as.numeric(n1$segment))
#    n1[,c('name','suffix'):=NULL]
#    n1$sample <- gsub('n1_','',n1$sample)
#    n2 <- result[grepl('n2_',name)]
#    n2[,suffix:=gsub('n2_','',name)]
#    n2_parsed <- rbindlist(lapply(unique(n2$name), parse))
#    n2 <- merge(n2, n2_parsed, by='name', all.x=T)
#    n2$segment <- round(as.numeric(n2$segment))
#    n2[,c('name','suffix'):=NULL]
#    n2$sample <- gsub('n2_','',n2$sample)
#    ascn <- merge(n1, n2, by=c('sample','segment','solution'))
#    setnames(ascn,c('value.x','value.y'),c('na','nb'))

    ## ploidy values for each solution
    pl <- data.table::dcast(name ~ solution, value.var='value', data=result[grepl('pl_',name),])
    pl <- cbind(sample=gsub('pl_','',pl$name), pl)
    pl[,name:=NULL]

    ## purity values for each solution
    pu <- data.table::dcast(name ~ solution, value.var='value', data=result[grepl('z_',name),])
    pu <- cbind(sample=gsub('z_','',pu$name), pu)
    pu[,name:=NULL]
    pu[,2:ncol(pu)] <- 1 / pu[,2:ncol(pu),with=F]

    ## 'allmatch' values for each segment, for each solution
    x <- data.table::dcast(name ~ solution, value.var='value', data=result[grepl('allmatch_',name),])
    x <- cbind(segment=round(as.numeric(gsub('allmatch_','',x$name))), x)
    x <- x[order(segment),]
    x[,name:=NULL]

    ## subset for UNIQUE solutions based only on allmatch, pl, and z values
    solutions <- names(x)[2:ncol(x)]

    if(length(solutions) > 1) {
        message('Reducing to unique solutions based on pu, pl, x values ...')
        current_pl <- pl[['sol0']]
        current_pu <- pu[['sol0']]
        current_x <- x[['sol0']]
        unique_solutions <- rep(as.logical(NA), length(solutions))
        names(unique_solutions) <- solutions
        unique_solutions['sol0'] <- T
        for(s in solutions[2:length(solutions)]) {
            next_pl <- pl[[s]]
            next_pu <- pu[[s]]
            next_x <- x[[s]]
            if(any(next_pl!=current_pl) | any(next_pu!=current_pu) | any(next_x!=current_x)) {
                unique_solutions[s] <- T  
            } else {
                unique_solutions[s] <- F    
            }
        }
        duplicate_solutions <- names(unique_solutions[unique_solutions==F])
        pl[,(duplicate_solutions):=NULL]
        pu[,(duplicate_solutions):=NULL]
        x[,(duplicate_solutions):=NULL]
        #ascn <- ascn[!solution %in% duplicate_solutions]
        message('Removed duplicate solutions: ',paste(duplicate_solutions,collapse=', '))
    }

    #ascn[,tcn:=na+nb]
    #ascn[na >= nb,mcn:=nb]
    #ascn[na < nb,mcn:=na]

    #pu <- cbind(variable='purity', pu)
    #pl <- cbind(variable='ploidy', pl)
    #x <- cbind(variable='allmatch', x)
    list(purity=pu, ploidy=pl, allmatch=x) #, ascn=ascn)
}


