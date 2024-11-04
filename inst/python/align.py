
import pandas as pd
import gurobipy as gb
from gurobipy import GRB

## dat should be a data.frame object from R with columns: "sample", "segment", "logR", "BAF", "GC", "mb"
def do_CNalign(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5.0, max_homdel_mb=50.0, delta=0.1, rho=0.85, both_alleles_must_align=1, gurobi_license='', epsilon=1e-4):

    # Create an environment with your WLS license
    with open(gurobi_license) as file:
        lines = [line.rstrip() for line in file]

    params = {
        "WLSACCESSID": lines[3].split('=')[1],
        "WLSSECRET": lines[4].split('=')[1],
        "LICENSEID": int(lines[5].split('=')[1]),
    }

    # Read input data into pandas DataFrame
    Samples = dat['sample'].unique()
    Segments = dat['segment'].unique()
    n_Samples = len(Samples)
    n_Segments = len(Segments)
    dat.set_index(['sample','segment'], inplace=True) ## set indices: sample, segment

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('ploidy range: ['+str(min_ploidy)+'-'+str(max_ploidy)+']')  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('Min length of any aligned segment (Mb): '+str(min_aligned_seg_mb))  
    print('Max hom-del length allowed (Mb): '+str(max_homdel_mb))  
    print('rho (max fraction of samples with matching segment): '+str(rho))  
    print('delta (max diff from integer CN): '+str(delta))  
    print('# samples in file: '+str(n_Samples))
    print('# segments in file: '+str(n_Segments))
    print('Both alleles must align: '+str(both_alleles_must_align))
    print('-------------------------------------')

    env = gb.Env(params=params)
    model = gb.Model(env=env)

    pl = {}
    z = {}  # z=1/pu
    for t in Samples:
        pl[t] = model.addVar(vtype=GRB.CONTINUOUS, name='pl_'+str(t), lb=min_ploidy, ub=max_ploidy)
        z[t] = model.addVar(vtype=GRB.CONTINUOUS, name='z_'+str(t), lb=1/max_purity, ub=1/min_purity)

    # Create auxilary continuous variables err s.t. err{t,s} = |y{t,s} - (x_{t,s}*m_{t}-b{t})|
    # Create binary variables match s.t. match{s}==1 if sum_s(err{t,s} <= delta) >= rho*num_samples
    err1 = {}
    err2 = {}
    n1 = {}
    n2 = {}
    n1int = {}
    n2int = {}
    neq0 = {}
    homdel_mb = {}
    n1bar = {}
    n2bar = {}
    large_enough = {}
    close_enough1 = {}
    close_enough2 = {}
    same_int1 = {}
    same_int2 = {}
    n1lwr = {}
    n1upr = {}
    n1cna = {}
    n2lwr = {}
    n2upr = {}
    n2cna = {}
    is_cna = {}
    match1 = {}
    match2 = {}
    match_both = {}
    match_both_with_cna = {}
    match_both_with_cna_and_large_enough = {}
    allmatch = {} 
    avg_ploidy = {}
    wt_copies1 = {} # segment-specific because sex-chromosomes can have different expected number of WT copies
    wt_copies2 = {} # segment-specific because sex-chromosomes can have different expected number of WT copies


    # variable for the average ploidy
    avg_ploidy = model.addVar(vtype=GRB.CONTINUOUS, name='avg_ploidy', lb=min_ploidy, ub=max_ploidy) # the average ploidy across all samples

    ## create segment-level variables
    for s in Segments:
        n1bar[s] = model.addVar(vtype=GRB.CONTINUOUS, name='n1bar_'+str(s))
        n2bar[s] = model.addVar(vtype=GRB.CONTINUOUS, name='n2bar_'+str(s))
        allmatch[s] = model.addVar(vtype=GRB.BINARY, name='allmatch_'+str(s))
        wt_copies1[s] = model.addVar(vtype=GRB.INTEGER, name='wt_copies1_'+str(s), lb=1, ub=6) # if both alleles have this number of copies, then don't count as CNA-ed
        wt_copies2[s] = model.addVar(vtype=GRB.INTEGER, name='wt_copies2_'+str(s), lb=1, ub=6) # if both alleles have this number of copies, then don't count as CNA-ed

        ## create segment,sample-level variables
        for t in Samples:
            n1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='n1_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # individual copy number cannot be negative
            n2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='n2_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # individual copy number cannot be negative
            n1int[t, s] = model.addVar(vtype=GRB.INTEGER, name='n1int_'+str(t)+','+str(s))
            n2int[t, s] = model.addVar(vtype=GRB.INTEGER, name='n2int_'+str(t)+','+str(s))
            err1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err1_'+str(t)+','+str(s))
            err2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err2_'+str(t)+','+str(s))
            match1[t, s] = model.addVar(vtype=GRB.BINARY, name='match1_'+str(t)+','+str(s))
            match2[t, s] = model.addVar(vtype=GRB.BINARY, name='match2_'+str(t)+','+str(s))
            match_both[t, s] = model.addVar(vtype=GRB.BINARY, name='match_both_'+str(t)+','+str(s))
            match_both_with_cna[t, s] = model.addVar(vtype=GRB.BINARY, name='match_both_with_cna_'+str(t)+','+str(s))
            match_both_with_cna_and_large_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='match_both_with_cna_and_large_enough'+str(t)+','+str(s))
            neq0[t, s] = model.addVar(vtype=GRB.BINARY, name='neq0_'+str(t)+','+str(s))
            n1lwr[t, s] = model.addVar(vtype="B",name='n1lwr_'+str(t)+','+str(s))
            n1upr[t, s] = model.addVar(vtype="B",name='n1upr_'+str(t)+','+str(s))
            n1cna[t, s] = model.addVar(vtype="B",name='n1cna_'+str(t)+','+str(s))
            n2lwr[t, s] = model.addVar(vtype="B",name='n2lwr_'+str(t)+','+str(s))
            n2upr[t, s] = model.addVar(vtype="B",name='n2upr_'+str(t)+','+str(s))
            n2cna[t, s] = model.addVar(vtype="B",name='n2cna_'+str(t)+','+str(s))
            is_cna[t, s] = model.addVar(vtype="B",name='is_cna_'+str(t)+','+str(s))
            large_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='large_enough_'+str(t)+','+str(s))
            close_enough1[t, s] = model.addVar(vtype=GRB.BINARY, name='close_enough1_'+str(t)+','+str(s))
            close_enough2[t, s] = model.addVar(vtype=GRB.BINARY, name='close_enough2_'+str(t)+','+str(s))
            same_int1[t, s] = model.addVar(vtype=GRB.BINARY, name='same_int1_'+str(t)+','+str(s))
            same_int2[t, s] = model.addVar(vtype=GRB.BINARY, name='same_int2_'+str(t)+','+str(s))

    
    ## get total homdel Mb in each sample
    for t in Samples:
        homdel_mb[t] = model.addVar(vtype=GRB.CONTINUOUS, name='homdel_mb_'+str(t), lb=0, ub=max_homdel_mb)

    ## the number of WT copies should be round(avg_ploidy/2)
    silence = model.addConstr(avg_ploidy == gb.quicksum(pl[t] for t in Samples)/n_Samples, name='c_pl_avg')

    ## segment,sample-level contraints
    for s in Segments:

        for t in Samples:
            ## calculate values
            r = dat.loc[t,s].logR
            b = dat.loc[t,s].BAF
            g = dat.loc[t,s].GC # germline copies
            c = 2**r
            c1 = 2**(r+1)

            ## check if segment is large enough
            l = dat.loc[t,s].mb
            silence = model.addGenConstrIndicator(large_enough[(t, s)], 1, l, GRB.GREATER_EQUAL, min_aligned_seg_mb, name='c_large_enough_'+str(t)+','+str(s))
 
            if b!=-9:
                ## logR+BAF are available so get allele-specific CN {n1,n2}
                silence = model.addConstr(n1[(t, s)] == -b*c*pl[t] + b*c1 - b*c1*z[t] + c*pl[t] - c1 + c1*z[t] + g - g*z[t] - 1 + z[t], name='c_n1_'+str(t)+','+str(s))
                silence = model.addConstr(n2[(t, s)] == b*c*pl[t] - b*c1 + b*c1*z[t] + 1 - z[t], name='c_n2_'+str(t)+','+str(s))

                ## check if copy1 differs from the nearest int by <=delta
                silence = model.addConstr(err1[(t, s)] >= n1int[(t,s)] - n1[(t,s)], name='c_err1_'+str(t)+','+str(s))
                silence = model.addConstr(err1[(t, s)] >= -n1int[(t,s)] + n1[(t,s)], name='c_err1neg_'+str(t)+','+str(s))
                silence = model.addGenConstrIndicator(close_enough1[(t, s)], 1, err1[(t, s)], GRB.LESS_EQUAL, delta, name='c_close_enough1_'+str(t)+','+str(s))
            
                ## check if copy2 differs from the nearest int by <=delta
                silence = model.addConstr(err2[(t, s)] >= n2int[(t,s)] - n2[(t,s)], name='c_err2_'+str(t)+','+str(s))
                silence = model.addConstr(err2[(t, s)] >= -n2int[(t,s)] + n2[(t,s)], name='c_err2neg_'+str(t)+','+str(s))
                silence = model.addGenConstrIndicator(close_enough2[(t, s)], 1, err2[(t, s)], GRB.LESS_EQUAL, delta, name='c_close_enough2_'+str(t)+','+str(s))
 
                ## constrain the rounded copy number to be within 0.5 of the *average* rounded copy number
                silence = model.addConstr(err1[(t, s)] >= n1bar[s] - n1int[(t, s)], name='c_err1_'+str(t)+','+str(s))
                silence = model.addConstr(err1[(t, s)] >= -n1bar[s] + n1int[(t, s)], name='c_err1neg_'+str(t)+','+str(s))     
                silence = model.addGenConstrIndicator(same_int1[(t, s)], 1, err1[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_same_int1_'+str(t)+','+str(s))
                silence = model.addGenConstrAnd(match1[(t,s)], [close_enough1[(t, s)], same_int1[(t,s)]], name='c_match1_'+str(t)+','+str(s)) 
                silence = model.addConstr(err2[(t, s)] >= n2bar[s] - n2int[(t, s)], name='c_err2_'+str(t)+','+str(s))
                silence = model.addConstr(err2[(t, s)] >= -n2bar[s] + n2int[(t, s)], name='c_err2neg_'+str(t)+','+str(s))     
                silence = model.addGenConstrIndicator(same_int2[(t, s)], 1, err2[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_same_int2_'+str(t)+','+str(s))
                silence = model.addGenConstrAnd(match2[(t,s)], [close_enough2[(t, s)], same_int2[(t,s)]], name='c_match2_'+str(t)+','+str(s)) 

                ## constraint for either n1 or n2 to have a CNA (not equal to the number of WT-copies)
                silence = model.addGenConstrOr(n1cna[(t, s)], [n1lwr[(t, s)], n1upr[(t, s)]], name='c_n1cna_'+str(t)+','+str(s)) 
                silence = model.addGenConstrOr(n2cna[(t, s)], [n2lwr[(t, s)], n2upr[(t, s)]], name='c_n2cna_'+str(t)+','+str(s)) 
                silence = model.addGenConstrOr(is_cna[(t, s)], [n1cna[(t, s)], n2cna[(t, s)]], name='c_is_cna_'+str(t)+','+str(s)) 

                ## constrain the rounded-n1,n2 variables
                silence = model.addConstr(n1int[(t,s)] <= n1[(t,s)] + 0.5 - epsilon, name='c_n1int_fwd_'+str(t)+','+str(s)) 
                silence = model.addConstr(n1int[(t,s)] >= n1[(t,s)] - 0.5, name='c_n1int_rev_'+str(t)+','+str(s)) 
                silence = model.addConstr(n2int[(t,s)] <= n2[(t,s)] + 0.5 - epsilon, name='c_n2int_fwd_'+str(t)+','+str(s)) 
                silence = model.addConstr(n2int[(t,s)] >= n2[(t,s)] - 0.5, name='c_n2int_rev_'+str(t)+','+str(s)) 

                ## constraint for CNAs in either allele
                silence = model.addConstr((n1upr[(t,s)]==1) >> (n1[(t, s)] >= wt_copies1[s]+0.5-epsilon))
                silence = model.addConstr((n1lwr[(t,s)]==1) >> (n1[(t, s)] <= wt_copies1[s]-0.5))
                silence = model.addConstr((n2upr[(t,s)]==1) >> (n2[(t, s)] >= wt_copies2[s]+0.5-epsilon))
                silence = model.addConstr((n2lwr[(t,s)]==1) >> (n2[(t, s)] <= wt_copies2[s]-0.5))

                ## check if it has TCN=0
                silence = model.addGenConstrIndicator(neq0[(t, s)], 1, n1[(t, s)] + n2[(t,s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_neq0_'+str(t)+','+str(s))

                ## check if each allele matches its nearest int
                if both_alleles_must_align==1:
                    silence = model.addGenConstrAnd(match_both[(t, s)], [match1[(t, s)], match2[(t, s)]], name='c_match_both_'+str(t)+','+str(s)) 
                else:
                    silence = model.addGenConstrOr(match_both[(t, s)], [match1[(t, s)], match2[(t, s)]], name='c_match_both_'+str(t)+','+str(s)) 

                if g==2:
                    silence = model.addConstr(wt_copies1[s] <= avg_ploidy/2 + 0.5 - epsilon, name='c_wt_copies1_fwd_'+str(s))
                    silence = model.addConstr(wt_copies1[s] >= avg_ploidy/2 - 0.5, name='c_wt_copies1_rev_'+str(s))
                    silence = model.addConstr(wt_copies2[s] <= avg_ploidy/2 + 0.5 - epsilon, name='c_wt_copies2_fwd_'+str(s))
                    silence = model.addConstr(wt_copies2[s] >= avg_ploidy/2 - 0.5, name='c_wt_copies2_rev_'+str(s))
                else:
                    silence = model.addConstr(wt_copies1[s] <= 1.5 - epsilon, name='c_wt_copies1_fwd_'+str(s))
                    silence = model.addConstr(wt_copies1[s] >= 0.5, name='c_wt_copies1_rev_'+str(s))
                    silence = model.addConstr(wt_copies2[s] == 0, name='c_wt_copies2'+str(s))

            else:
                # only logR is available, so get total CN {n1, NA} s.t. n1=total number of copies
                silence = model.addConstr(n1[(t, s)] == z[t]*c*2 + c*pl[t] - 2*c - g*z[t] + g, name='c_n1_'+str(t)+','+str(s))
                silence = model.addConstr(n2[(t, s)] == 0, name='c_n2_'+str(t)+','+str(s))

                ## check if copy1 differs from the nearest int by <=delta
                silence = model.addConstr(err1[(t, s)] >= n1int[(t,s)] - n1[(t,s)], name='c_err1_'+str(t)+','+str(s))
                silence = model.addConstr(err1[(t, s)] >= -n1int[(t,s)] + n1[(t,s)], name='c_err1neg_'+str(t)+','+str(s))
                silence = model.addGenConstrIndicator(close_enough1[(t, s)], 1, err1[(t, s)], GRB.LESS_EQUAL, delta, name='c_close_enough1_'+str(t)+','+str(s))
            
                ## constrain the rounded copy number to be within 0.5 of the *average* rounded copy number
                silence = model.addConstr(err1[(t, s)] >= n1bar[s] - n1int[(t, s)], name='c_err1_'+str(t)+','+str(s))
                silence = model.addConstr(err1[(t, s)] >= -n1bar[s] + n1int[(t, s)], name='c_err1neg_'+str(t)+','+str(s))     
                silence = model.addGenConstrIndicator(same_int1[(t, s)], 1, err1[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_same_int1_'+str(t)+','+str(s))
                silence = model.addGenConstrAnd(match1[(t,s)], [close_enough1[(t, s)], same_int1[(t,s)]], name='c_match1_'+str(t)+','+str(s)) 

                ## constraint for n1 (only) to have a CNA (not equal to the number of WT-copies)
                silence = model.addGenConstrOr(n1cna[(t, s)], [n1lwr[(t, s)], n1upr[(t, s)]], name='c_n1cna_'+str(t)+','+str(s)) 
                silence = model.addConstr(is_cna[(t, s)] == n1cna[(t, s)], name='c_is_cna_'+str(t)+','+str(s)) 

                ## constrain the rounded-n1 variable
                silence = model.addConstr(n1int[(t,s)] <= n1[(t,s)] + 0.5 - epsilon, name='c_n1int_fwd_'+str(t)+','+str(s)) 
                silence = model.addConstr(n1int[(t,s)] >= n1[(t,s)] - 0.5, name='c_n1int_rev_'+str(t)+','+str(s)) 
                silence = model.addConstr(n2int[(t,s)] == 0, name='c_n2int_'+str(t)+','+str(s)) 

                ## constraint for CNAs in either allele
                silence = model.addConstr((n1upr[(t,s)]==1) >> (n1[(t, s)] >= wt_copies1[s]+0.5-epsilon))
                silence = model.addConstr((n1lwr[(t,s)]==1) >> (n1[(t, s)] <= wt_copies1[s]-0.5))

                ## check if the only allele (n1) matches its nearest int
                silence = model.addConstr(match_both[(t, s)] == match1[(t, s)], name='c_match_both_'+str(t)+','+str(s)) 

                if g==2:
                    silence = model.addConstr(wt_copies1[s] <= avg_ploidy + 0.5 - epsilon, name='c_wt_copies1_fwd_'+str(s))
                    silence = model.addConstr(wt_copies1[s] >= avg_ploidy - 0.5, name='c_wt_copies1_rev_'+str(s))
                else:
                    silence = model.addConstr(wt_copies1[s] <= 1.5 - epsilon, name='c_wt_copies1_fwd_'+str(s))
                    silence = model.addConstr(wt_copies1[s] >= 0.5, name='c_wt_copies1_rev_'+str(s))
 
            ## segment,sample-level constraint that both alleles satisfy conditions to match other samples (and also have CNAs, so we don't just return diploid segments)
            silence = model.addGenConstrAnd(match_both_with_cna[(t, s)], [match_both[(t, s)], is_cna[(t, s)]], name='c_match_both_with_cna_'+str(t)+','+str(s)) 
            silence = model.addGenConstrAnd(match_both_with_cna_and_large_enough[(t, s)], [match_both_with_cna[(t, s)], large_enough[(t,s)]], name='c_match_both_with_cna_and_large_enough_'+str(t)+','+str(s)) 



    for t in Samples:
        ## sample-level constraints
        silence = model.addConstr(homdel_mb[t] == gb.quicksum(dat.loc[(t,s)].mb * neq0[(t, s)] for s in Segments), name='c_homdel_mb_'+str(t))

    ## the hard part is to check that the nearest integer for each copy n1int,n2int is the SAME integer for at least rho% of samples
    ## NB: this should work for both BAF available (n2!=0) and TCN-only (n2=0)
    for s in Segments:

        ## get the average of the rounded n1,n2
        silence = model.addConstr(n1bar[s] == gb.quicksum(n1int[(t, s)] for t in Samples)/n_Samples, name='c_n1bar_'+str(s))
        silence = model.addConstr(n2bar[s] == gb.quicksum(n2int[(t, s)] for t in Samples)/n_Samples, name='c_n2bar_'+str(s))
        
        # create a variable for the objective function, which is a binary variable for whether 
        # a segment was a 'match' in at least rho% of samples.
        model.addGenConstrIndicator(allmatch[s], 1, gb.quicksum(match_both_with_cna_and_large_enough[(t, s)] for t in Samples), GRB.GREATER_EQUAL, rho*n_Samples, name='c_allmatch_'+str(s))

    model.setObjective(gb.quicksum(allmatch[s] for s in Segments), GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return model



## non-allele specific version.
## dat should be a data.frame object from R with columns: "sample", "segment", "logR", "GC"
def CNalign_tcn(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_homdels=0, max_homdels=3, t_diff=0.05, t_err=0.1, rho=0.85, gurobi_license='', epsilon=1e-4, assume_wgd=False):

    # Create an environment with your WLS license
    with open(gurobi_license) as file:
        lines = [line.rstrip() for line in file]

    params = {
        "WLSACCESSID": lines[3].split('=')[1],
        "WLSSECRET": lines[4].split('=')[1],
        "LICENSEID": int(lines[5].split('=')[1]),
    }

    # Read input data into pandas DataFrame
    Samples = dat['sample'].unique()
    Segments = dat['segment'].unique()
    n_Samples = len(Samples)
    n_Segments = len(Segments)
    dat.set_index(['sample','segment'], inplace=True) ## set indices: sample, segment

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('ploidy range: ['+str(min_ploidy)+'-'+str(max_ploidy)+']')  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('Total hom-dels allowed: ['+str(min_homdels)+'-'+str(max_homdels)+']')  
    print('rho (max fraction of samples with matching segment): '+str(rho))  
    print('t_diff (max diff from integer CN): '+str(t_diff))  
    print('t_err (max diff from integer CN): '+str(t_err))  
    print('# samples in file: '+str(n_Samples))
    print('# segments in file: '+str(n_Segments))
    print('-------------------------------------')

    env = gb.Env(params=params)
    model = gb.Model(env=env)

    pl = {}
    z = {}  # z=1/pu
    for t in Samples:
        pl[t] = model.addVar(vtype=GRB.CONTINUOUS, name='pl_'+str(t), lb=min_ploidy, ub=max_ploidy)
        z[t] = model.addVar(vtype=GRB.CONTINUOUS, name='z_'+str(t), lb=1/max_purity, ub=1/min_purity)

    # Create auxilary continuous variables err s.t. err{t,s} = |y{t,s} - (x_{t,s}*m_{t}-b{t})|
    # Create binary variables match s.t. match{s}==1 if sum_s(err{t,s} <= threshold) == rho*num_samples
    err = {}
    diff = {}
    n = {}
    n_int = {}
    n_eq0 = {}
    n_homdel = {}
    n_bar = {}
    n_lwr = {}
    n_upr = {}
    n_cna = {}
    is_cna = {}
    match_diff = {}
    match_err = {}
    match = {}
    match_with_cna = {}
    allmatch = {}

    ## create segment-level variables
    for s in Segments:
        n_bar[s] = model.addVar(vtype=GRB.CONTINUOUS, name='n_bar_'+str(s))
        allmatch[s] = model.addVar(vtype=GRB.BINARY, name='allmatch_'+str(s))

        ## create segment,sample-level variables
        for t in Samples:
            n[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='n_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # total copy number cannot be negative
            n_int[t, s] = model.addVar(vtype=GRB.INTEGER, name='n_int_'+str(t)+','+str(s))
            err[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err_'+str(t)+','+str(s))
            diff[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='diff_'+str(t)+','+str(s))
            match_diff[t, s] = model.addVar(vtype=GRB.BINARY, name='match_diff_'+str(t)+','+str(s))
            match_err[t, s] = model.addVar(vtype=GRB.BINARY, name='match_err_'+str(t)+','+str(s))
            match[t, s] = model.addVar(vtype=GRB.BINARY, name='match_'+str(t)+','+str(s))
            match_with_cna[t, s] = model.addVar(vtype=GRB.BINARY, name='match_with_cna'+str(t)+','+str(s))
            n_eq0[t, s] = model.addVar(vtype=GRB.BINARY, name='n_eq0_'+str(t)+','+str(s))
            n_lwr[t, s] = model.addVar(vtype="B",name='n_lwr_'+str(t)+','+str(s))
            n_upr[t, s] = model.addVar(vtype="B",name='n_upr_'+str(t)+','+str(s))
            n_cna[t, s] = model.addVar(vtype="B",name='n_cna_'+str(t)+','+str(s))
            is_cna[t, s] = model.addVar(vtype="B",name='is_cna_'+str(t)+','+str(s))

    ## get number of homdel segments in each sample
    for t in Samples:
        n_homdel[t] = model.addVar(vtype=GRB.INTEGER, name='n_homdel_'+str(t), lb=min_homdels, ub=max_homdels)

    ## segment,sample-level contraints
    for s in Segments:
        for t in Samples:
            ## calculate values
            r = dat.loc[t,s].logR
            g = dat.loc[t,s].GC # germline copies
            c = 2**r
            c1 = 2**(r+1)

            # only logR is available, so get total CN {n, NA} s.t. n=total number of copies
            silence = model.addConstr(n[(t, s)] == z[t]*c*2 + c*pl[t] - 2*c - g*z[t] + g, name='c_n_'+str(t)+','+str(s))
            silence = model.addConstr(n[(t, s)] - 0.5 + epsilon <= n_int[(t,s)], name='c_n_int_1_'+str(t)+','+str(s))
            silence = model.addConstr(n_int[(t,s)] - 0.5 <= n[(t,s)], name='c_n_int_2_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(n_eq0[(t, s)], 1, n[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_n_eq0_'+str(t)+','+str(s))

    for t in Samples:
        ## sample-level constraints
        silence = model.addConstr(n_homdel[t] == gb.quicksum(n_eq0[(t, s)] for s in Segments), name='c_n_homdel_'+str(t))

    for s in Segments:
        ## segment-level constraints
        silence = model.addConstr(n_bar[s] == gb.quicksum(n[(t, s)] for t in Samples), name='c_n_bar_'+str(s))
    
        for t in Samples:
            ## segment,sample-level constraints
            ## for allele1
            silence = model.addConstr(diff[(t, s)] >= n_int[(t,s)] - n[(t, s)], name='c_diff_'+str(t)+','+str(s))
            silence = model.addConstr(diff[(t, s)] >= -n_int[(t,s)] + n[(t, s)], name='c_diffneg_'+str(t)+','+str(s))     
            silence = model.addConstr(err[(t, s)] >= n_bar[s]/n_Samples - n[(t, s)], name='c_err_'+str(t)+','+str(s))
            silence = model.addConstr(err[(t, s)] >= -n_bar[s]/n_Samples + n[(t, s)], name='c_errneg_'+str(t)+','+str(s))     
            silence = model.addGenConstrIndicator(match_diff[(t, s)], 1, diff[(t, s)], GRB.LESS_EQUAL, t_diff, name='c_match_diff_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(match_err[(t, s)], 1, err[(t, s)], GRB.LESS_EQUAL, t_err, name='c_match_err_'+str(t)+','+str(s))
            silence = model.addConstr(match[(t, s)] == match_diff[(t,s)]*match_err[(t, s)], name='c_match_'+str(t)+','+str(s)) ## needs to be ==
            ## constraint for TCN to have a CNA (not equal 2, or 4 for WGD cases)
            if assume_wgd is False:
                silence = model.addGenConstrIndicator(n_upr[(t,s)], 1, n[(t, s)], GRB.GREATER_EQUAL, 2.5) 
                silence = model.addGenConstrIndicator(n_lwr[(t,s)], 1, n[(t, s)], GRB.LESS_EQUAL, 1.5-epsilon) 
            else:
                silence = model.addGenConstrIndicator(n_upr[(t,s)], 1, n[(t, s)], GRB.GREATER_EQUAL, 4.5) 
                silence = model.addGenConstrIndicator(n_lwr[(t,s)], 1, n[(t, s)], GRB.LESS_EQUAL, 3.5-epsilon) 
          
            silence = model.addGenConstrOr(n_cna[(t, s)], [n_lwr[(t, s)], n_upr[(t, s)]], name='c_n_cna_'+str(t)+','+str(s)) 
            ## segment,sample-level constraint that the segment satisfy conditions to match other samples (and also have CNAs, so we don't just return diploid/tetraploid segments)
            silence = model.addGenConstrAnd(match_with_cna[(t, s)], [match[(t, s)], n_cna[(t, s)]], name='c_match_with_cna'+str(t)+','+str(s)) 

    # create a variable for the objective function, which is a binary variable for whether 
    # a segment was a 'match' in at least rho% of samples.
    for s in Segments:
        model.addGenConstrIndicator(allmatch[s], 1, gb.quicksum(match_with_cna[(t, s)] for t in Samples), GRB.GREATER_EQUAL, rho*n_Samples, name='c_allmatch_'+str(s))

    model.setObjective(gb.quicksum(allmatch[s] for s in Segments), GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return model


