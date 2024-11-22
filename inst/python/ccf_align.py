
import pandas as pd
import gurobipy as gb
from gurobipy import GRB

## dat should be a data.frame object from R with columns: "sample", "variant", "vaf", "gc"
def do_CCFalign(dat, min_purity=0.05, max_purity=0.95, max_tcn=6, min_ccf_is_clonal=0.9, min_frac_clonal_is_truncal=1.0, gurobi_license=''):

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
    Variants = dat['variant'].unique()
    n_Samples = len(Samples)
    n_Variants = len(Variants)
    dat.set_index(['sample','variant'], inplace=True) ## set indices: sample, variant

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('max TCN: '+str(max_tcn))
    print('min CCF considered clonal: '+str(min_ccf_is_clonal))
    print('min fraction of clonal samples for truncal variant: '+str(min_frac_clonal_is_truncal))
    print('# samples in file: '+str(n_Samples))
    print('# variants in file: '+str(n_Variants))
    print('-------------------------------------')

    env = gb.Env(params=params)
    model = gb.Model(env=env)

    pu = {}
    tcn = {}
    mcn = {}
    ccf = {}
    LHS1 = {}
    LHS2 = {}
    RHS = {}
    clonal = {} # variant s clonal in sample t
    fraction_clonal = {} # fraction of tumors where sample is clonal
    truncal = {} # variant is clonal in every sample

    for t in Samples:
        pu[t] = model.addVar(vtype=GRB.CONTINUOUS, name='pu_'+str(t), lb=min_purity, ub=max_purity)

    ## create variant-level variables
    for s in Variants:
        fraction_clonal[s] = model.addVar(vtype=GRB.CONTINUOUS, name='fraction_clonal_'+str(s), lb=0, ub=1) 
        truncal[s] = model.addVar(vtype=GRB.BINARY, name='truncal_'+str(s)) # is fraction_clonal >= min_fraction_clonal?

        ## create variant,sample-level variables
        for t in Samples:
            tcn[t, s] = model.addVar(vtype=GRB.INTEGER, name='tcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            mcn[t, s] = model.addVar(vtype=GRB.INTEGER, name='mcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # MCN
            ccf[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='ccf_'+str(t)+','+str(s), lb=0, ub=1) # CCF
            LHS1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS1_'+str(t)+','+str(s), lb=0, ub=1) # eq. LHS
            LHS2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS2_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. LHS
            RHS[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='RHS_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. RHS
            clonal[t, s] = model.addVar(vtype=GRB.BINARY, name='clonal_'+str(t)+','+str(s)) # is variant clonal in sample t?

    ## variant-level contraints
    for s in Variants:
        silence = model.addConstr(fraction_clonal[s] == gb.quicksum(clonal[(t, s)] for t in Samples)/n_Samples, name='c_fraction_clonal_'+str(s))
        silence = model.addGenConstrIndicator(truncal[s], 1, fraction_clonal[s], GRB.GREATER_EQUAL, min_frac_clonal_is_truncal, name='c_truncal_'+str(s))
 
        ## variant,sample-level contraints
        for t in Samples:
            vaf = dat.loc[t,s].vaf # vaf
            g = dat.loc[t,s].gc # germline copies
            silence = model.addConstr(mcn[(t,s)] <= tcn[(t,s)], name='c_mcn_leq_tcn_'+str(t)+','+str(s))
            silence = model.addConstr(LHS1[(t,s)] == pu[t]*ccf[(t,s)], name='c_LHS1_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == LHS1[(t,s)]*mcn[(t,s)], name='c_LHS2_'+str(t)+','+str(s))
            silence = model.addConstr(RHS[(t,s)] == (pu[t]*vaf*tcn[(t,s)] + (1-pu[t])*g), name='c_RHS_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == RHS[(t,s)], name='c_LHS2_EQ_RHS_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(clonal[(t, s)], 1, ccf[(t,s)], GRB.GREATER_EQUAL, min_ccf_is_clonal, name='c_clonal_'+str(t)+','+str(s))

    model.setObjective(gb.quicksum(truncal[s] for s in Variants), GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return model


## dat should be a data.frame object from R with columns: "sample", "variant", "vaf", "gc"
def do_CCFalign2(dat, min_purity=0.05, max_purity=0.95, max_tcn=6, gurobi_license=''):

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
    Variants = dat['variant'].unique()
    n_Samples = len(Samples)
    n_Variants = len(Variants)
    dat.set_index(['sample','variant'], inplace=True) ## set indices: sample, variant

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('max TCN: '+str(max_tcn))
    #print('min CCF considered clonal: '+str(min_ccf_is_clonal))
    #print('min fraction of clonal samples for truncal variant: '+str(min_frac_clonal_is_truncal))
    print('# samples in file: '+str(n_Samples))
    print('# variants in file: '+str(n_Variants))
    print('-------------------------------------')

    env = gb.Env(params=params)
    model = gb.Model(env=env)

    pu = {}
    tcn = {}
    mcn = {}
    ccf = {}
    LHS1 = {}
    LHS2 = {}
    RHS = {}
    avgCCFmut = {} 

    for t in Samples:
        pu[t] = model.addVar(vtype=GRB.CONTINUOUS, name='pu_'+str(t), lb=min_purity, ub=max_purity)

    ## create variant-level variables
    for s in Variants:
        avgCCFmut[s] = model.addVar(vtype=GRB.CONTINUOUS, name='avgCCFmut_'+str(s), lb=0, ub=1) # average CCF/mut across samples

        ## create variant,sample-level variables
        for t in Samples:
            tcn[t, s] = model.addVar(vtype=GRB.INTEGER, name='tcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            mcn[t, s] = model.addVar(vtype=GRB.INTEGER, name='mcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # MCN
            ccf[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='ccf_'+str(t)+','+str(s), lb=0, ub=1) # CCF
            LHS1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS1_'+str(t)+','+str(s), lb=0, ub=1) # eq. LHS
            LHS2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS2_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. LHS
            RHS[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='RHS_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. RHS

    ## variant-level contraints
    for s in Variants:
        silence = model.addConstr(avgCCFmut[s] == gb.quicksum(ccf[(t,s)] for t in Samples)/n_Samples, name='c_avgCCFmut_'+str(s))
        
        ## variant,sample-level contraints
        for t in Samples:
            vaf = dat.loc[t,s].vaf # vaf
            g = dat.loc[t,s].gc # germline copies
            silence = model.addConstr(mcn[(t,s)] <= tcn[(t,s)], name='c_mcn_leq_tcn_'+str(t)+','+str(s))
            silence = model.addConstr(LHS1[(t,s)] == pu[t]*ccf[(t,s)], name='c_LHS1_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == LHS1[(t,s)]*mcn[(t,s)], name='c_LHS2_'+str(t)+','+str(s))
            silence = model.addConstr(RHS[(t,s)] == (pu[t]*vaf*tcn[(t,s)] + (1-pu[t])*g), name='c_RHS_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == RHS[(t,s)], name='c_LHS2_EQ_RHS_'+str(t)+','+str(s))

    ## our objective is to maximize the overall average CCF across all mutations and samples
    model.setObjective(gb.quicksum(avgCCFmut[s] for s in Variants) / n_Variants, GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return model


