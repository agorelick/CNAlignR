##' run_CNAlign
##'
##' Determine purity and ploidy values for multiple tumor samples with shared ancestry. This function uses the GuRoBi solver to determine purity/ploidy values for each sample that will *maximize* the number of segments with the same (allele-specific) integer copy numbers in at least rho% of samples. 
##'
##' @export
run_CNAlign <- function(dat, gurobi_license, py_script=NA, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5.0, max_homdel_mb=100.0, 
                    delta_tcn_to_int=0.2, delta_tcn_to_avg=0.1, delta_tcnavg_to_int=0.1, delta_mcn_to_int=0.2, delta_mcn_to_avg=0.1, delta_mcnavg_to_int=0.1, 
                    rho=1.0, timeout=3*60, min_cna_segments_per_sample=1, mcn_weight=0.5, obj2_clonalonly=F, sol_count=10, max_tcn_avg_int=100, tcn_only=F, n_unique_tcn_in_sol=1) {
                
    require(lubridate)
    require(reticulate)
    if(is.na(py_script)) py_script <- system.file('CNAlign/CNAlign/core.py',package='CNAlignR')
	source_python(py_script)

    # recode NA BAFs to -9
	dat$BAF[is.na(dat$BAF)] <- -9
	
	start_time <- now()
	message('Started CNAlign at ',as.character(start_time))

    if(tcn_only) {
        message('Running in TCN-only mode.')
        dat$BAF <- -9
        mcn_weight <- 0
        delta_mcn_to_int <- 0.5
        delta_mcn_to_avg <- 0.5 
        delta_mcnavg_to_int <- 0.5
    }

    df <- CNAlign(dat, gurobi_license, min_ploidy, max_ploidy, min_purity, max_purity, 
                  min_aligned_seg_mb, max_homdel_mb, 
                  delta_tcn_to_int, delta_tcn_to_avg, delta_tcnavg_to_int, 
                  delta_mcn_to_int, delta_mcn_to_avg, delta_mcnavg_to_int, 
                  mcn_weight, rho, timeout, min_cna_segments_per_sample, obj2_clonalonly, sol_count, max_tcn_avg_int, n_unique_tcn_in_sol)

	end_time <- now()
	run_date <- as.character(format(end_time,format='%Y-%m-%d %H:%M'))
    message('Ended CNAlign at ',as.character(end_time))
    sec_elapsed <- round(as.numeric(end_time) - as.numeric(start_time))
    message('Time elapsed: ',sec_elapsed,'s')

    params <- list(min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, min_aligned_seg_mb=min_aligned_seg_mb, max_homdel_mb=max_homdel_mb, delta_tcn_to_int=delta_tcn_to_int, delta_tcn_to_avg=delta_tcn_to_avg, delta_tcnavg_to_int=delta_tcnavg_to_int, delta_mcn_to_int=delta_mcn_to_int, delta_mcn_to_avg=delta_mcn_to_avg, delta_mcnavg_to_int=delta_mcnavg_to_int, rho=rho, gurobi_license=license, timeout=timeout, min_cna_segments_per_sample=min_cna_segments_per_sample, mcn_weight=mcn_weight, obj2_clonalonly=obj2_clonalonly, sol_count=sol_count, max_tcn_avg_int=max_tcn_avg_int)

    # return the model object and a list of all the parameters
    out <- list(df=df, params=params)
}


parse_brackets <- function(text) { 
    require(stringr)
    before <- str_extract(text, ".*?(?=\\[)")
    inside <- str_extract(text, "(?<=\\[).*?(?=\\])")
    after <- str_replace(text, ".*\\] ?", "")
    list(text=text, before=before, inside=inside, after=after)
}



get_solutions <- function(m, num_solutions=NA, ncpus=1) { 
    require(parallel)
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
        vars <- rbindlist(lapply(values_for_solution$name, parse_brackets))
        values_for_solution <- cbind(values_for_solution, vars)
        values_for_solution[!is.na(before), category:=before]
        values_for_solution[is.na(before), category:=text]
        values_for_solution[!is.na(inside), subject:=inside]
        message('Solution ',solnumber,' extracted.')
        values_for_solution[,c('text','before','inside','after'):=NULL]
        values_for_solution
    }
    if(is.na(num_solutions)) {
        solutions <- 0:m$SolCount
    } else {
        solutions <- 0:(num_solutions-1)   
    }

    message('Extracting values for ',length(solutions),' solutions ...')
    result_list <- mclapply(solutions, get_values_for_solution, m, mc.cores=ncpus)
    result <- rbindlist(result_list, fill=T)
    result$solution <- paste0('sol',result$solution)
    solution_labels <- unique(result$solution)
    result$solution <- factor(result$solution, levels=solution_labels)
    result[category=='z', category:='pu']
    result[category=='pu', value:=1/value]
    result <- data.table::dcast(name + category + subject ~ solution, value.var='value', data=result)

    # extract purity solutions
    pu <- result[category=='pu',c('subject',solution_labels),with=F]
    pl <- result[category=='pl',c('subject',solution_labels),with=F]
    allmatch <- result[category=='allmatch',c('subject',solution_labels),with=F]

    # subset for unique solutions
    if(length(solution_labels) > 1) {
        message('Reducing to unique solutions based on pu[sample], pl[sample], allmatch[segment] values ...')
        current_pl <- pl[['sol0']]
        current_pu <- pu[['sol0']]
        current_allmatch <- allmatch[['sol0']]
        unique_solutions <- rep(as.logical(NA), length(solutions))
        names(unique_solutions) <- solution_labels
        unique_solutions['sol0'] <- T
        for(s in solution_labels[2:length(solutions)]) {
            next_pl <- pl[[s]]
            next_pu <- pu[[s]]
            next_allmatch <- allmatch[[s]]
            if(any(next_pl!=current_pl) | any(next_pu!=current_pu) | any(next_allmatch!=current_allmatch)) {
                unique_solutions[s] <- T  
            } else {
                unique_solutions[s] <- F    
            }
        }
        duplicate_solutions <- names(unique_solutions[unique_solutions==F])
        pl[,(duplicate_solutions):=NULL]
        pu[,(duplicate_solutions):=NULL]
        allmatch[,(duplicate_solutions):=NULL]
    }

    sol0_highlight <- as.integer(gsub('seg','',allmatch[allmatch$sol0==1]$subject))
    sol0_pu <- pu[,c('subject','sol0'),with=F]
    sol0_pu$variable <- 'pu'
    sol0_pl <- pl[,c('subject','sol0'),with=F]
    sol0_pl$variable <- 'pl'
    sol0_fits <- rbind(sol0_pu, sol0_pl)
    sol0_fits <- data.table::dcast(subject ~ variable, value.var='sol0', data=sol0_fits)        
    names(sol0_fits)[1] <- 'sample'
    list(sol0_fits=sol0_fits, sol0_highlight=sol0_highlight, pu=pu, pl=pl, allmatch=allmatch)
}



