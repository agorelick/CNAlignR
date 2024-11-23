
import pandas as pd
import gurobipy as gb
from gurobipy import GRB

## dat should be a data.frame object from R with columns: "sample", "variant", "vaf", "gc"
def do_CCFalign(dat, min_purity=0.05, max_purity=0.95, max_tcn=6, rho=0.8, delta=0.1, gurobi_license=''):

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
    print('rho: '+str(rho))
    print('delta: '+str(delta))
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
    tcn_int = {}
    mcn_int = {}
    tcn_err = {}
    mcn_err = {}
    tcn_near_int = {}
    mcn_near_int = {}
    LHS1 = {}
    LHS2 = {}
    RHS = {}
    avgccf = {} # avg CCF for each variant (across the samples)
    avgmcn = {} # avg MCN for each variant (across the samples)
    mcn_diff_from_avg = {} # for each variant this is its difference from the average MCN across all samples. When this is 0 then all samples have the same MCN
    mcn_close_enough = {} 
    mcn_near_int_and_close_enough = {}
    mcn_frac_match = {}
    mcn_matched = {}
    avgtcn = {} # avg TCN for each variant (across the samples)
    tcn_diff_from_avg = {} # for each variant this is its difference from the average TCN across all samples. When this is 0 then all samples have the same TCN
    tcn_close_enough = {} 
    tcn_frac_match = {}
    tcn_matched = {}
    tcn_near_int_and_close_enough = {}
    ascn_matched = {}

    for t in Samples:
        pu[t] = model.addVar(vtype=GRB.CONTINUOUS, name='pu_'+str(t), lb=min_purity, ub=max_purity)

    ## create variant-level variables
    for s in Variants:
        avgccf[s] = model.addVar(vtype=GRB.CONTINUOUS, name='avgccf_'+str(s), lb=0, ub=1) 
        avgmcn[s] = model.addVar(vtype=GRB.CONTINUOUS, name='avgmcn_'+str(s), lb=1, ub=max_tcn) 
        avgtcn[s] = model.addVar(vtype=GRB.CONTINUOUS, name='avgtcn_'+str(s), lb=1.7, ub=max_tcn) 
        mcn_frac_match[s] = model.addVar(vtype=GRB.CONTINUOUS, name='mcn_frac_match_'+str(s),lb=0,ub=1) 
        tcn_frac_match[s] = model.addVar(vtype=GRB.CONTINUOUS, name='tcn_frac_match_'+str(s),lb=0,ub=1) 
        mcn_matched[s] = model.addVar(vtype=GRB.BINARY, name='mcn_matched_'+str(s)) 
        tcn_matched[s] = model.addVar(vtype=GRB.BINARY, name='tcn_matched_'+str(s)) 
        ascn_matched[s] = model.addVar(vtype=GRB.BINARY, name='ascn_matched_'+str(s)) 

        ## create variant,sample-level variables
        for t in Samples:
            tcn[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='tcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            tcn_int[t, s] = model.addVar(vtype=GRB.INTEGER, name='tcn_int_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            tcn_err[t, s] = model.addVar(vtype=GRB.INTEGER, name='tcn_err_'+str(t)+','+str(s), lb=0, ub=0.5) # TCN
            tcn_near_int[t, s] = model.addVar(vtype=GRB.BINARY, name='tcn_near_int_'+str(t)+','+str(s)) # TCN
            tcn_diff_from_avg[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='tcn_diff_from_avg_'+str(t)+','+str(s)) 
            tcn_close_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='tcn_close_enough_'+str(t)+','+str(s)) 
            tcn_near_int_and_close_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='tcn_near_int_and_close_enough_'+str(t)+','+str(s)) 

            mcn[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='mcn_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            mcn_int[t, s] = model.addVar(vtype=GRB.INTEGER, name='mcn_int_'+str(t)+','+str(s), lb=1, ub=max_tcn) # TCN
            mcn_err[t, s] = model.addVar(vtype=GRB.INTEGER, name='mcn_err_'+str(t)+','+str(s), lb=0, ub=0.5) # TCN
            mcn_near_int[t, s] = model.addVar(vtype=GRB.BINARY, name='mcn_near_int_'+str(t)+','+str(s)) # TCN
            mcn_diff_from_avg[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='mcn_diff_from_avg_'+str(t)+','+str(s)) 
            mcn_close_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='mcn_close_enough_'+str(t)+','+str(s)) 
            mcn_near_int_and_close_enough[t, s] = model.addVar(vtype=GRB.BINARY, name='mcn_near_int_and_close_enough_'+str(t)+','+str(s)) 

            ccf[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='ccf_'+str(t)+','+str(s), lb=0, ub=1) # CCF
            LHS1[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS1_'+str(t)+','+str(s), lb=0, ub=1) # eq. LHS
            LHS2[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='LHS2_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. LHS
            RHS[t, s] = model.addVar(vtype=GRB.CONTINUOUS, name='RHS_'+str(t)+','+str(s), lb=0, ub=max_tcn) # eq. RHS

    ## variant-level contraints
    for s in Variants:

        silence = model.addConstr(avgccf[s] == gb.quicksum(ccf[(t, s)] for t in Samples)/n_Samples, name='c_avgccf_'+str(s))
        silence = model.addConstr(avgtcn[s] == gb.quicksum(tcn[(t, s)] for t in Samples)/n_Samples, name='c_avgtcn_'+str(s))
        silence = model.addConstr(avgmcn[s] == gb.quicksum(mcn[(t, s)] for t in Samples)/n_Samples, name='c_avgmcn_'+str(s))

        ## to consider a variant 'matched' we require that a min fraction (rho) of samples have MCN near integer value and close to the overall avg MCN.
        silence = model.addConstr(tcn_frac_match[s] == gb.quicksum(tcn_near_int_and_close_enough[(t,s)] for t in Samples)/n_Samples, name='c_tcn_frac_match_'+str(s))
        silence = model.addConstr(mcn_frac_match[s] == gb.quicksum(mcn_near_int_and_close_enough[(t,s)] for t in Samples)/n_Samples, name='c_mcn_frac_match_'+str(s))
        silence = model.addGenConstrIndicator(mcn_matched[s], 1, mcn_frac_match[s], GRB.GREATER_EQUAL, rho, name='c_mcn_matched_'+str(s))
        silence = model.addGenConstrIndicator(tcn_matched[s], 1, tcn_frac_match[s], GRB.GREATER_EQUAL, rho, name='c_tcn_matched_'+str(s))
        silence = model.addGenConstrAnd(ascn_matched[s], [mcn_matched[s], tcn_matched[s]], name='c_ascn_matched_'+str(s))

        ## variant,sample-level contraints
        for t in Samples:
            vaf = dat.loc[t,s].vaf # vaf
            g = dat.loc[t,s].gc # germline copies

            silence = model.addGenConstrIndicator(mcn_near_int[(t, s)], 1, mcn_err[(t, s)], GRB.LESS_EQUAL, 0.5, name='c_mcn_near_int_'+str(t)+','+str(s))
            silence = model.addConstr(mcn[(t,s)] <= tcn[(t,s)], name='c_mcn_leq_tcn_'+str(t)+','+str(s))
            silence = model.addConstr(mcn_diff_from_avg[t, s] >= mcn[(t,s)] - avgmcn[s], name='c_mcn_diff_from_avg_fwd_'+str(t)+','+str(s))
            silence = model.addConstr(mcn_diff_from_avg[t, s] >= -mcn[(t,s)] + avgmcn[s], name='c_mcn_diff_from_avg_rev_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(mcn_close_enough[t, s], 1, mcn_diff_from_avg[t, s], GRB.LESS_EQUAL, delta, name='c_mcn_close_enough_'+str(t)+','+str(s))
            silence = model.addGenConstrAnd(mcn_near_int_and_close_enough[t, s], [mcn_near_int[(t,s)], mcn_close_enough[(t, s)]], name='c_mcn_near_int_and_close_enough_'+str(t)+','+str(s))

            silence = model.addGenConstrIndicator(tcn_near_int[(t, s)], 1, tcn_err[(t, s)], GRB.LESS_EQUAL, 0.5, name='c_tcn_near_int_'+str(t)+','+str(s))
            silence = model.addConstr(tcn_diff_from_avg[t, s] >= tcn[(t,s)] - avgtcn[s], name='c_tcn_diff_from_avg_fwd_'+str(t)+','+str(s))
            silence = model.addConstr(tcn_diff_from_avg[t, s] >= -tcn[(t,s)] + avgtcn[s], name='c_tcn_diff_from_avg_rev_'+str(t)+','+str(s))
            silence = model.addGenConstrIndicator(tcn_close_enough[t, s], 1, tcn_diff_from_avg[t, s], GRB.LESS_EQUAL, delta, name='c_tcn_close_enough_'+str(t)+','+str(s))
            silence = model.addGenConstrAnd(tcn_near_int_and_close_enough[t, s], [tcn_near_int[(t,s)], tcn_close_enough[(t, s)]], name='c_tcn_near_int_and_close_enough_'+str(t)+','+str(s))

            silence = model.addConstr(LHS1[(t,s)] == pu[t]*ccf[(t,s)], name='c_LHS1_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == LHS1[(t,s)]*mcn[(t,s)], name='c_LHS2_'+str(t)+','+str(s))
            silence = model.addConstr(RHS[(t,s)] == (pu[t]*vaf*tcn[(t,s)] + (1-pu[t])*g), name='c_RHS_'+str(t)+','+str(s))
            silence = model.addConstr(LHS2[(t,s)] == RHS[(t,s)], name='c_LHS2_EQ_RHS_'+str(t)+','+str(s))
            silence = model.addConstr(mcn_err[(t,s)] >= mcn[(t,s)] - mcn_int[(t,s)], name='c_mcn_err_fwd_'+str(s))
            silence = model.addConstr(mcn_err[(t,s)] >= -mcn[(t,s)] + mcn_int[(t,s)], name='c_mcn_err_rev_'+str(s))
            silence = model.addConstr(tcn_err[(t,s)] >= tcn[(t,s)] - tcn_int[(t,s)], name='c_tcn_err_fwd_'+str(s))
            silence = model.addConstr(tcn_err[(t,s)] >= -tcn[(t,s)] + tcn_int[(t,s)], name='c_tcn_err_rev_'+str(s))


    model.setObjective(gb.quicksum(avgccf[s]*ascn_matched[s] for s in Variants)/n_Variants, GRB.MAXIMIZE)
    model.update()
    model.optimize()
    return model



