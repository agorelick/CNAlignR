
import pandas as pd
import gurobipy as gb
from gurobipy import GRB

## dat should be a data.frame object from R with columns: "sample", "segment", "logR", "BAF", "GC"
def CNalign(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_homdels=0, max_homdels=3, t_diff=0.05, t_err=0.1, rho=0.85, both_alleles_must_align=1, gurobi_license='', epsilon=1e-4, assume_wgd=False):

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
    print('# hom-dels allowed: ['+str(min_homdels)+'-'+str(max_homdels)+']')  
    print('rho (max fraction of samples with matching segment): '+str(rho))  
    print('t_diff (max diff from integer CN): '+str(t_diff))  
    print('t_err (max diff from integer CN): '+str(t_err))  
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
    # Create binary variables match s.t. match{s}==1 if sum_s(err{t,s} <= threshold) == rho*num_samples
    err1 = {}
    err2 = {}
    err3 = {}
    diff1 = {}
    diff2 = {}
    diff3 = {}
    n1 = {}
    n2 = {}
    n1int = {}
    n2int = {}
    neq0 = {}
    n_homdel = {}
    n1bar = {}
    n2bar = {}
    tcn = {}
    n1lwr = {}
    n1upr = {}
    n1cna = {}
    n2lwr = {}
    n2upr = {}
    n2cna = {}
    is_cna = {}
    match_diff1 = {}
    match_err1 = {}
    match_diff2 = {}
    match_err2 = {}
    match1 = {}
    match2 = {}
    match_both = {}
    match_both_with_cna = {}
    allmatch = {}

    ## create segment-level variables
    for s in Segments:
        n1bar[s] = model.addVar(vtype=GRB.CONTINUOUS, name='n1bar_'+str(s))
        n2bar[s] = model.addVar(vtype=GRB.CONTINUOUS, name='n2bar_'+str(s))
        allmatch[s] = model.addVar(vtype=GRB.BINARY, name='allmatch_'+str(s))

        ## create segment,sample-level variables
        for t in Samples:
            n1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='n1_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # individual copy number cannot be negative
            n2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='n2_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # individual copy number cannot be negative
            tcn[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='tcn_'+str(t)+','+str(s), lb=-0.5+epsilon, ub=1e6) # total copy number cannot be negative
            n1int[t, s] = model.addVar(vtype=GRB.INTEGER, name='n1int_'+str(t)+','+str(s))
            n2int[t, s] = model.addVar(vtype=GRB.INTEGER, name='n2int_'+str(t)+','+str(s))
            err1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err1_'+str(t)+','+str(s))
            err2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err2_'+str(t)+','+str(s))
            err3[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='err3_'+str(t)+','+str(s))
            diff1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='diff1_'+str(t)+','+str(s))
            diff2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='diff2_'+str(t)+','+str(s))
            diff3[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='diff3_'+str(t)+','+str(s))
            match_diff1[t, s] = model.addVar(vtype=GRB.BINARY, name='match_diff1_'+str(t)+','+str(s))
            match_err1[t, s] = model.addVar(vtype=GRB.BINARY, name='match_err1_'+str(t)+','+str(s))
            match_diff2[t, s] = model.addVar(vtype=GRB.BINARY, name='match_diff2_'+str(t)+','+str(s))
            match_err2[t, s] = model.addVar(vtype=GRB.BINARY, name='match_err2_'+str(t)+','+str(s))
            match1[t, s] = model.addVar(vtype=GRB.BINARY, name='match1_'+str(t)+','+str(s))
            match2[t, s] = model.addVar(vtype=GRB.BINARY, name='match2_'+str(t)+','+str(s))
            match_both[t, s] = model.addVar(vtype=GRB.BINARY, name='match_both_'+str(t)+','+str(s))
            match_both_with_cna[t, s] = model.addVar(vtype=GRB.BINARY, name='match_both_with_cna'+str(t)+','+str(s))
            neq0[t, s] = model.addVar(vtype=GRB.BINARY, name='neq0_'+str(t)+','+str(s))
            n1lwr[t, s] = model.addVar(vtype="B",name='n1lwr_'+str(t)+','+str(s))
            n1upr[t, s] = model.addVar(vtype="B",name='n1upr_'+str(t)+','+str(s))
            n1cna[t, s] = model.addVar(vtype="B",name='n1cna_'+str(t)+','+str(s))
            n2lwr[t, s] = model.addVar(vtype="B",name='n2lwr_'+str(t)+','+str(s))
            n2upr[t, s] = model.addVar(vtype="B",name='n2upr_'+str(t)+','+str(s))
            n2cna[t, s] = model.addVar(vtype="B",name='n2cna_'+str(t)+','+str(s))
            is_cna[t, s] = model.addVar(vtype="B",name='is_cna_'+str(t)+','+str(s))

    ## get number of homdel segments in each sample
    for t in Samples:
        n_homdel[t] = model.addVar(vtype=GRB.INTEGER, name='n_homdel_'+str(t), lb=min_homdels, ub=max_homdels)

    ## segment,sample-level contraints
    for s in Segments:
        for t in Samples:
            ## calculate values
            r = dat.loc[t,s].logR
            b = dat.loc[t,s].BAF
            g = dat.loc[t,s].GC # germline copies
            c = 2**r
            c1 = 2**(r+1)

            if pd.isna(b)==False:
                ## logR+BAF are available so get allele-specific CN {n1,n2}
                silence = model.addConstr(n1[(t, s)] == -b*c*pl[t] + b*c1 - b*c1*z[t] + c*pl[t] - c1 + c1*z[t] + g - g*z[t] - 1 + z[t], name='c_n1_'+str(t)+','+str(s))
                silence = model.addConstr(n1[(t, s)] - 0.5 + epsilon <= n1int[(t,s)], name='c_n1int_1_'+str(t)+','+str(s))
                silence = model.addConstr(n1int[(t,s)] - 0.5 <= n1[(t,s)], name='c_n1int_2_'+str(t)+','+str(s))
                silence = model.addConstr(n2[(t, s)] == b*c*pl[t] - b*c1 + b*c1*z[t] + 1 - z[t], name='c_n2_'+str(t)+','+str(s))
                silence = model.addConstr(n2[(t, s)] - 0.5 + epsilon <= n2int[(t,s)], name='c_n2int_1_'+str(t)+','+str(s))
                silence = model.addConstr(n2int[(t,s)] - 0.5 <= n2[(t,s)], name='c_n2int_2_'+str(t)+','+str(s))
                silence = model.addConstr(tcn[(t,s)] == n1[(t,s)] + n2[(t,s)], name='c_tcn_'+str(t)+','+str(s))
                silence = model.addGenConstrIndicator(neq0[(t, s)], 1, n1[(t, s)] + n2[(t,s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_neq0_'+str(t)+','+str(s))
            else:
                # only logR is available, so get total CN {n1, NA} s.t. n1=total number of copies
                silence = model.addConstr(n1[(t, s)] == z[t]*c*2 + c*pl[t] - 2*c - g*z[t] + g, name='c_n1_'+str(t)+','+str(s))
                silence = model.addConstr(n1[(t, s)] - 0.5 + epsilon <= n1int[(t,s)], name='c_n1int_1_'+str(t)+','+str(s))
                silence = model.addConstr(n1int[(t,s)] - 0.5 <= n1[(t,s)], name='c_n1int_2_'+str(t)+','+str(s))
                silence = model.addConstr(n2[(t, s)] == 0, name='c_n2_'+str(t)+','+str(s)) # we are forcing n1+n2=TCN and n2=0
                silence = model.addGenConstrIndicator(neq0[(t, s)], 1, n1[(t, s)] + n2[(t,s)], GRB.LESS_EQUAL, 0.5-epsilon, name='c_neq0_'+str(t)+','+str(s))
                silence = model.addConstr(tcn[(t,s)] == n1[(t,s)] + n2[(t,s)], name='c_tcn_'+str(t)+','+str(s))

    for t in Samples:
        ## sample-level constraints
        silence = model.addConstr(n_homdel[t] == gb.quicksum(neq0[(t, s)] for s in Segments), name='c_n_homdel_'+str(t))

    for s in Segments:
        ## segment-level constraints
        silence = model.addConstr(n1bar[s] == gb.quicksum(n1[(t, s)] for t in Samples), name='c_n1bar_'+str(s))
        silence = model.addConstr(n2bar[s] == gb.quicksum(n2[(t, s)] for t in Samples), name='c_n2bar_'+str(s))
    
        for t in Samples:
            ## segment,sample-level constraints
            ## for allele1
            silence = model.addConstr(diff1[(t, s)] >= n1int[(t,s)] - n1[(t, s)], name='c_diff1_'+str(t)+','+str(s))
            silence = model.addConstr(diff1[(t, s)] >= -n1int[(t,s)] + n1[(t, s)], name='c_diff1neg_'+str(t)+','+str(s))     
            silence = model.addConstr(err1[(t, s)] >= n1bar[s]/n_Samples - n1[(t, s)], name='c_err1_'+str(t)+','+str(s))
            silence = model.addConstr(err1[(t, s)] >= -n1bar[s]/n_Samples + n1[(t, s)], name='c_err1neg_'+str(t)+','+str(s))     
            silence = model.addGenConstrIndicator(match_diff1[(t, s)], 1, diff1[(t, s)], GRB.LESS_EQUAL, t_diff, name='c_match_diff1_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(match_err1[(t, s)], 1, err1[(t, s)], GRB.LESS_EQUAL, t_err, name='c_match_err1_'+str(t)+','+str(s))
            silence = model.addConstr(match1[(t, s)] == match_diff1[(t,s)]*match_err1[(t, s)], name='c_match1_'+str(t)+','+str(s)) ## needs to be ==
            ## for allele2
            silence = model.addConstr(diff2[(t, s)] >= n2int[(t,s)] - n2[(t, s)], name='c_diff2_'+str(t)+','+str(s))
            silence = model.addConstr(diff2[(t, s)] >= -n2int[(t,s)] + n2[(t, s)], name='c_diff2neg_'+str(t)+','+str(s))     
            silence = model.addConstr(err2[(t, s)] >= n2bar[s]/n_Samples - n2[(t, s)], name='c_err2_'+str(t)+','+str(s))
            silence = model.addConstr(err2[(t, s)] >= -n2bar[s]/n_Samples + n2[(t, s)], name='c_err2neg_'+str(t)+','+str(s))     
            silence = model.addGenConstrIndicator(match_diff2[(t, s)], 1, diff2[(t, s)], GRB.LESS_EQUAL, t_diff, name='c_match_diff2_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(match_err2[(t, s)], 1, err2[(t, s)], GRB.LESS_EQUAL, t_err, name='c_match_err2_'+str(t)+','+str(s))
            silence = model.addConstr(match2[(t, s)] == match_diff2[(t,s)]*match_err2[(t, s)], name='c_match2_'+str(t)+','+str(s)) ## needs to be ==
            ## constraint for either n1 or n2 to have a CNA (not equal 1)
            if assume_wgd is False:
                silence = model.addGenConstrIndicator(n1upr[(t,s)], 1, n1[(t, s)], GRB.GREATER_EQUAL, 1.5) 
                silence = model.addGenConstrIndicator(n1lwr[(t,s)], 1, n1[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon) 
                silence = model.addGenConstrIndicator(n2upr[(t,s)], 1, n2[(t, s)], GRB.GREATER_EQUAL, 1.5) 
                silence = model.addGenConstrIndicator(n2lwr[(t,s)], 1, n2[(t, s)], GRB.LESS_EQUAL, 0.5-epsilon) 
            else:
                silence = model.addGenConstrIndicator(n1upr[(t,s)], 1, n1[(t, s)], GRB.GREATER_EQUAL, 2.5) 
                silence = model.addGenConstrIndicator(n1lwr[(t,s)], 1, n1[(t, s)], GRB.LESS_EQUAL, 1.5-epsilon) 
                silence = model.addGenConstrIndicator(n2upr[(t,s)], 1, n2[(t, s)], GRB.GREATER_EQUAL, 2.5) 
                silence = model.addGenConstrIndicator(n2lwr[(t,s)], 1, n2[(t, s)], GRB.LESS_EQUAL, 1.5-epsilon) 
          
            silence = model.addGenConstrOr(n1cna[(t, s)], [n1lwr[(t, s)], n1upr[(t, s)]], name='c_n1cna_'+str(t)+','+str(s)) 
            silence = model.addGenConstrOr(n2cna[(t, s)], [n2lwr[(t, s)], n2upr[(t, s)]], name='c_n2cna_'+str(t)+','+str(s)) 
            silence = model.addGenConstrOr(is_cna[(t, s)], [n1cna[(t, s)], n2cna[(t, s)]], name='c_is_cna_'+str(t)+','+str(s)) 
            if both_alleles_must_align==1:
                silence = model.addGenConstrAnd(match_both[(t, s)], [match1[(t, s)], match2[(t, s)]], name='c_match_both_'+str(t)+','+str(s)) 
            else:
                silence = model.addGenConstrOr(match_both[(t, s)], [match1[(t, s)], match2[(t, s)]], name='c_match_both_'+str(t)+','+str(s)) 
            ## segment,sample-level constraint that both alleles satisfy conditions to match other samples (and also have CNAs, so we don't just return diploid segments)
            silence = model.addGenConstrAnd(match_both_with_cna[(t, s)], [match_both[(t, s)], is_cna[(t, s)]], name='c_match_both_with_cna'+str(t)+','+str(s)) 

    # create a variable for the objective function, which is a binary variable for whether 
    # a segment was a 'match' in at least rho% of samples.
    for s in Segments:
        model.addGenConstrIndicator(allmatch[s], 1, gb.quicksum(match_both_with_cna[(t, s)] for t in Samples), GRB.GREATER_EQUAL, rho*n_Samples, name='c_allmatch_'+str(s))

    model.setObjective(gb.quicksum(allmatch[s] for s in Segments), GRB.MAXIMIZE)
    model.update()
    model.write('1939.lp')
    model.optimize()
    return model


