##' CCFalign
##' @export
CCFalign <- function(dat, min_purity=0.05, max_purity=0.95, max_tcn=6, min_ccf_is_clonal=0.9, min_frac_clonal_is_truncal=1.0, gurobi_license='') {
    require(reticulate)
    require(lubridate)
    if(is.na(py_script)) py_script <- system.file("python", "align.py", package = "CNalign")

    ## dat should be a data.frame object from R with columns: "sample", "variant", "vaf", "gc"
	source_python(py_script)
	start_time <- now()
	message('Started CNalign at ',as.character(start_time))
	m <- do_CCFalign(dat, min_purity=min_purity, max_purity=max_purity, max_tcn=max_tcn, min_ccf_is_clonal=min_ccf_is_clonal, min_frac_clonal_is_truncal=min_frac_clonal_is_truncal, gurobi_license=gurobi_license) 
	end_time <- now()
	run_date <- as.character(format(end_time,format='%Y-%m-%d %H:%M'))
    message('Ended CNalign at ',as.character(end_time))
    sec_elapsed <- round(as.numeric(end_time) - as.numeric(start_time))
    message('Time elapsed: ',sec_elapsed,'s')

    params <- list(min_purity=min_purity, max_purity=max_purity, max_tcn=max_tcn, min_ccf_is_clonal=min_ccf_is_clonal, min_frac_clonal_is_truncal=min_frac_clonal_is_truncal, gurobi_license=gurobi_license) 

    # return the model object and a list of all the parameters
    out <- list(m=m, params=params)
}


##' CNalign
##'
##' Determine purity and ploidy values for multiple tumor samples with shared ancestry. This function uses the GuRoBi solver to determine purity/ploidy values for each sample that will *maximize* the number of segments with the same (allele-specific) integer copy numbers in at least rho% of samples. 
##'
##' @export
CNalign <- function(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5, max_homdel_mb=20, delta=0.1, rho=0.85, both_alleles_must_align=1, epsilon=1e-4, tcn_only=F, gurobi_license='~/gurobi.lic', py_script=NA, aligned_includes_wt=0) {
    require(reticulate)
    require(lubridate)
    if(is.na(py_script)) py_script <- system.file("python", "align.py", package = "CNalign")

    ## given input matrices of aligned sample/segment:
    ## 1. logR (required, no NAs)
    ## 2. BAF (optional, can be fully/partially NAs)
    ## Determine purity/ploidy values for each sample such that there is maximum number of segments with the same copy number values for >= r% of the samples

	dat$BAF[is.na(dat$BAF)] <- -9
	if(tcn_only==T) {
		message('TCN-only alignment')
		dat$BAF <- -9 # remove all BAFs to force TCN-only alignment
		both_alleles_must_align <- 0
	}
	
	source_python(py_script)
	start_time <- now()
	message('Started CNalign at ',as.character(start_time))
	m <- do_CNalign(dat, min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_aligned_seg_mb=min_aligned_seg_mb, max_homdel_mb=max_homdel_mb, delta=delta, rho=rho, both_alleles_must_align=both_alleles_must_align, gurobi_license=gurobi_license, epsilon=epsilon, aligned_includes_wt=aligned_includes_wt)
	end_time <- now()
	run_date <- as.character(format(end_time,format='%Y-%m-%d %H:%M'))
    message('Ended CNalign at ',as.character(end_time))
    sec_elapsed <- round(as.numeric(end_time) - as.numeric(start_time))
    message('Time elapsed: ',sec_elapsed,'s')

    params <- list(min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_aligned_seg_mb=min_aligned_seg_mb, max_homdel_mb=max_homdel_mb, delta=delta, rho=rho, both_alleles_must_align=both_alleles_must_align, gurobi_license=gurobi_license, epsilon=epsilon, tcn_only=tcn_only, run_date=run_date, aligned_includes_wt=aligned_includes_wt)

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
	result <- data.table::dcast(name ~ solution, value.var='value', data=result)

	# extract purity solutions
	pu <- d2m(result[grepl('z_',name),])
	pu <- 1 / pu
	rownames(pu) <- gsub('z_','pu_',rownames(pu))
	pu <- cbind(name=rownames(pu), as.data.table(pu))

	# ploidy solutions
	pl <- result[grepl('pl_',name),]

	# OTHER VARIANTS
	other <- result[!name %in% c(pu$name,pl$name),] 	
	result <- rbind(pu, pl, other)

	# return everything
	result
}


get_pu_pl_solutions <- function(result) {
	## ploidy values for each solution
    pl <- result[grepl('pl_',name),]
    pl <- cbind(sample=gsub('pl_','',pl$name), pl)
    pl[,name:=NULL]

    ## purity values for each solution
    pu <- result[grepl('z_',name),]
    pu <- cbind(sample=gsub('z_','',pu$name), pu)
    pu[,name:=NULL]
    pu[,2:ncol(pu)] <- 1 / pu[,2:ncol(pu),with=F]

    ## 'allmatch' values for each segment, for each solution
    x <- result[grepl('allmatch_seg',name),]
    x <- cbind(segment=round(as.numeric(gsub('allmatch_seg','',x$name))), x)
    x <- x[order(segment),]
    x[,name:=NULL]

    ## subset for UNIQUE solutions based only on allmatch, pl, and z values
    solutions <- names(x)[2:ncol(x)]

    if(length(solutions) > 1) {
        message('Reducing to unique solutions based on pu, pl, x values ...')
        current_pl <- pl[['0']]
        current_pu <- pu[['0']]
        current_x <- x[['0']]
        unique_solutions <- rep(as.logical(NA), length(solutions))
        names(unique_solutions) <- solutions
        unique_solutions['0'] <- T
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
    list(purity=pu, ploidy=pl, allmatch=x)
}


