

# =============================================================================
# allele-specific, multi-objective:
# first, try to maximize the number of segments that align with CNAs in at at least one allele
# then, for those segments, try to get as close as possible to integer copy number in all samples
# =============================================================================

import pandas as pd
import gurobipy as gb
from gurobipy import GRB
import time
import numpy as np

## dat should be a data.frame object from R with columns: "sample", "segment", "logR", "BAF", "GC", "mb"
def do_CNalign(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5.0, max_homdel_mb=50.0, 
                delta_tcn_to_int=0.2, delta_tcn_to_avg=0.1, delta_tcnavg_to_int=0.05, delta_mcn_to_int=0.2, delta_mcn_to_avg=0.1, delta_mcnavg_to_int=0.05, mcn_weight=0.5,
                rho=1.0, gurobi_license='', normal_baseline=2, aligned_includes_wt=0, timeout=10*60, min_cna_segments_per_sample=3, obj2_clonalonly=0):

    # Create an environment with your WLS license
    with open(gurobi_license) as file:
        lines = [line.rstrip() for line in file]

    params = {
        "WLSACCESSID": lines[3].split('=')[1],
        "WLSSECRET": lines[4].split('=')[1],
        "LICENSEID": int(lines[5].split('=')[1]),
    }

    class StagnationCallback:
        def __init__(self, max_stagnation_seconds):
            self.max_stagnation = max_stagnation_seconds
            self.obj_count = 0
            self.last_improve_time = time.time()
            self.best_obj = float('inf')
            self.stagnated = False

        def __call__(self, model, where):
            now = time.time()

            if where == GRB.Callback.MULTIOBJ:
                # Starting a new multi-objective level
                self.obj_count = model.cbGet(GRB.Callback.MULTIOBJ_OBJCNT)
                self.best_obj = float('inf')
                self.last_improve_time = now
                self.stagnated = False
                print(f"[Callback] Starting objective {self.obj_count}")

            elif where == GRB.Callback.MIP:
                try:
                    current_obj = model.cbGet(GRB.Callback.MIP_OBJBST)
                    if current_obj < GRB.INFINITY:
                        if abs(current_obj - self.best_obj) > 1e-5:
                            self.best_obj = current_obj
                            self.last_improve_time = now
                        elif now - self.last_improve_time > self.max_stagnation and not self.stagnated:
                            print(f"[Callback] No improvement for {self.max_stagnation}s on objective {self.obj_count}. Moving to next.")
                            self.stagnated = True
                            model.cbStopOneMultiObj(self.obj_count)
                except Exception as e:
                    print(f"[Callback] Exception in MIP callback: {e}") 

    # create custom callback
    callback = StagnationCallback(max_stagnation_seconds=timeout)

    # Read input data into pandas DataFrame
    Samples = dat['sample'].unique()
    Segments = dat['segment'].unique()
    n_Samples = len(Samples)
    n_Segments = len(Segments)
    dat.set_index(['sample','segment'], inplace=True) ## set indices: sample, segment

    if np.isin(normal_baseline,[2,4,6])==False:
        return('Invalid baseline!')

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('ploidy range: ['+str(min_ploidy)+'-'+str(max_ploidy)+']')  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('Min length of any aligned segment (Mb): '+str(min_aligned_seg_mb))  
    print('Max hom-del length allowed (Mb): '+str(max_homdel_mb))  
    print('rho (min fraction of samples with matching segment): '+str(rho))  
    #print('delta (max diff from integer CN): '+str(delta))  
    print('# samples in file: '+str(n_Samples))
    print('# segments in file: '+str(n_Segments))
    #print('Wild-type copy numbers included in alignment: '+str(aligned_includes_wt))
    print('-------------------------------------')

    env = gb.Env(params=params)
    
    
    model = gb.Model(env=env)
            
    n1 = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='n1') 
    #n2 = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='n2') 

    tcn = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='tcn')    
    tcn_avg = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='tcn_avg')
    tcn_int = model.addVars(Samples, Segments, vtype=GRB.INTEGER, name='tcn_int')
    tcn_int_err = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='tcn_int_err', lb=0)
    tcn_spread = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='tcn_spread', lb=0)
    tcn_avg_int = model.addVars(Samples, Segments, vtype=GRB.INTEGER, name='tcn_avg_int', lb=0)
    tcn_avg_int_err = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='tcn_avg_int_err', lb=0)

    tcn_close_to_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_close_to_int')
    tcn_close_to_avg = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_close_to_avg')
    tcn_avg_close_to_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_avg_close_to_int')
    tcn_match = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_match')
    tcn_match_and_avg_at_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_match_and_avg_at_int')
    tcn_gain = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_gain')
    tcn_loss = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_loss')
    tcn_cna = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='tcn_cna')
    tcn_error_clonal = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='tcn_error_clonal')


    mcn = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='mcn')    
    mcn_avg = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='mcn_avg')
    mcn_int = model.addVars(Samples, Segments, vtype=GRB.INTEGER, name='mcn_int')
    mcn_int_err = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='mcn_int_err', lb=0)
    mcn_spread = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='mcn_spread', lb=0)
    mcn_avg_int = model.addVars(Samples, Segments, vtype=GRB.INTEGER, name='mcn_avg_int', lb=0)
    mcn_avg_int_err = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='mcn_avg_int_err', lb=0)

    mcn_close_to_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_close_to_int')
    mcn_close_to_avg = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_close_to_avg')
    mcn_avg_close_to_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_avg_close_to_int')
    mcn_match = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_match')    
    mcn_match_and_avg_at_int = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_match_and_avg_at_int')
    mcn_gain = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_gain')
    mcn_loss = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_loss')
    mcn_cna = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='mcn_cna')
    mcn_error_clonal = model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='mcn_error_clonal')
    
    
    # additional Sample+Segment-level variables
    match_both = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_both')
    match_both_and_large_enough = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_both_and_large_enough')        
    match_both_and_large_enough_and_cna = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_both_and_large_enough_and_cna')        
    is_homdel = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='is_homdel')
    is_cna = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='is_cna')
    large_enough = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='large_enough')
    
    # additional Segment-level variables
    allmatch = model.addVars(Segments, vtype=GRB.BINARY, name='allmatch')
    
    # additional Sample-level variables
    pl = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='pl', lb=min_ploidy, ub=max_ploidy)
    z = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='z', lb=1/max_purity, ub=1/min_purity)
    homdel_mb = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='homdel_mb', lb=0, ub=max_homdel_mb)
    n_cna_segments_per_sample = model.addVars(Samples, vtype=GRB.INTEGER, name='n_cna_segments_per_sample', lb=min_cna_segments_per_sample, ub=n_Segments)

    # variable for the average ploidy
    avg_ploidy = model.addVar(vtype=GRB.CONTINUOUS, name='avg_ploidy', lb=min_ploidy, ub=max_ploidy) # the average ploidy across all samples
        
    ## the number of WT copies should be round(avg_ploidy/2)
    model.addConstr(avg_ploidy == gb.quicksum(pl[t] for t in Samples)/n_Samples, name='c_pl_avg')

    # objective variables
    n_clonal = model.addVar(vtype=GRB.INTEGER, lb=0, ub=n_Segments, name='n_clonal')




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
            model.addGenConstrIndicator(large_enough[(t, s)], 1, l, GRB.GREATER_EQUAL, min_aligned_seg_mb, name='c_large_enough_'+str(t)+','+str(s))

            if b!=-9:
                ## logR+BAF are available so get n1 and n2
                model.addConstr(n1[t, s] == -b*c*pl[t] + b*c1 - b*c1*z[t] + c*pl[t] - c1 + c1*z[t] + g - g*z[t] - 1 + z[t])
                model.addConstr(mcn[t, s] == b*c*pl[t] - b*c1 + b*c1*z[t] + 1 - z[t])
                tcn_wt_copies = g
                mcn_wt_copies = g-1
            else:
                # only logR is available, so get total CN {n1, NA} s.t. n1=total number of copies
                model.addConstr(n1[(t, s)] == z[t]*c*2 + c*pl[t] - 2*c - g*z[t] + g, name='c_n1_'+str(t)+','+str(s))
                model.addConstr(mcn[(t, s)] == 0, name='c_n2_'+str(t)+','+str(s))
                tcn_wt_copies = g
                mcn_wt_copies = 0

            model.addConstr(tcn[t, s] == n1[t,s] + mcn[t,s])
            
            # =============================================================================
            # TCN
            # =============================================================================
            
            # is TCN close to its nearest integer
            model.addConstr(tcn_int[t,s] <= tcn[t,s] + 0.5) 
            model.addConstr(tcn_int[t,s] >= tcn[t,s] - 0.5) 
            model.addConstr(tcn_int_err[t,s] >= tcn_int[t,s] - tcn[t,s])
            model.addConstr(tcn_int_err[t,s] >= -tcn_int[t,s] + tcn[t,s])
            model.addGenConstrIndicator(tcn_close_to_int[t,s], 1, tcn_int_err[t,s], GRB.LESS_EQUAL, delta_tcn_to_int)
            
            # is TCN close to the TCNavg (not too spread out)
            #model.addConstr(tcn_avg[t,s] == gb.quicksum(tcn_int[t,s] for t in Samples)/n_Samples)
            #model.addConstr(tcn_spread[t,s] >= tcn_avg[t,s] - tcn_int[t,s])
            #model.addConstr(tcn_spread[t,s] >= -tcn_avg[t,s] + tcn_int[t,s])     
            model.addConstr(tcn_avg[t,s] == gb.quicksum(tcn[t,s] for t in Samples)/n_Samples)
            model.addConstr(tcn_spread[t,s] >= tcn_avg[t,s] - tcn[t,s])
            model.addConstr(tcn_spread[t,s] >= -tcn_avg[t,s] + tcn[t,s])     
            model.addGenConstrIndicator(tcn_close_to_avg[t,s], 1, tcn_spread[t,s], GRB.LESS_EQUAL, delta_tcn_to_avg)
            
            # is TCNavg close to its nearest integer
            model.addConstr(tcn_avg_int[t,s] <= tcn_avg[t,s] + 0.5) 
            model.addConstr(tcn_avg_int[t,s] >= tcn_avg[t,s] - 0.5) 
            model.addConstr(tcn_avg_int_err[t,s] >= tcn_avg[t,s] - tcn_avg_int[t,s])
            model.addConstr(tcn_avg_int_err[t,s] >= -tcn_avg[t,s] + tcn_avg_int[t,s])                
            model.addGenConstrIndicator(tcn_avg_close_to_int[t,s], 1, tcn_avg_int_err[t,s], GRB.LESS_EQUAL, delta_tcnavg_to_int)
            
            
            ## match if both close enough and same int as the rounded average
            model.addGenConstrAnd(tcn_match[t,s], [tcn_close_to_int[t,s], tcn_close_to_avg[t,s]]) 
            model.addGenConstrAnd(tcn_match_and_avg_at_int[t,s], [tcn_match[t,s], tcn_avg_close_to_int[t,s]]) 

            ## constraint for TCN-based CNA
            model.addConstr((tcn_gain[t,s]==1) >> (tcn_int[t,s] >= tcn_wt_copies + 1))
            model.addConstr((tcn_loss[t,s]==1) >> (tcn_int[t,s] <= tcn_wt_copies - 1))
            model.addGenConstrOr(tcn_cna[t,s], [tcn_gain[t,s], tcn_loss[t,s]])


            # =============================================================================
            # MCN
            # =============================================================================

            # is MCN close to its nearest integer
            model.addConstr(mcn_int[t,s] <= mcn[t,s] + 0.5) 
            model.addConstr(mcn_int[t,s] >= mcn[t,s] - 0.5) 
            model.addConstr(mcn_int_err[t,s] >= mcn_int[t,s] - mcn[t,s])
            model.addConstr(mcn_int_err[t,s] >= -mcn_int[t,s] + mcn[t,s])
            model.addGenConstrIndicator(mcn_close_to_int[t,s], 1, mcn_int_err[t,s], GRB.LESS_EQUAL, delta_mcn_to_int)
            
            # is MCN close to the MCNavg (not too spread out)
            #model.addConstr(mcn_avg[t,s] == gb.quicksum(mcn_int[t,s] for t in Samples)/n_Samples)
            #model.addConstr(mcn_spread[t,s] >= mcn_avg[t,s] - mcn_int[t,s])
            #model.addConstr(mcn_spread[t,s] >= -mcn_avg[t,s] + mcn_int[t,s])     
            model.addConstr(mcn_avg[t,s] == gb.quicksum(mcn[t,s] for t in Samples)/n_Samples)
            model.addConstr(mcn_spread[t,s] >= mcn_avg[t,s] - mcn[t,s])
            model.addConstr(mcn_spread[t,s] >= -mcn_avg[t,s] + mcn[t,s])     
            model.addGenConstrIndicator(mcn_close_to_avg[t,s], 1, mcn_spread[t,s], GRB.LESS_EQUAL, delta_mcn_to_avg)
            
            # is MCNavg close to its nearest integer
            model.addConstr(mcn_avg_int[t,s] <= mcn_avg[t,s] + 0.5) 
            model.addConstr(mcn_avg_int[t,s] >= mcn_avg[t,s] - 0.5) 
            model.addConstr(mcn_avg_int_err[t,s] >= mcn_avg[t,s] - mcn_avg_int[t,s])
            model.addConstr(mcn_avg_int_err[t,s] >= -mcn_avg[t,s] + mcn_avg_int[t,s])                
            model.addGenConstrIndicator(mcn_avg_close_to_int[t,s], 1, mcn_avg_int_err[t,s], GRB.LESS_EQUAL, delta_mcnavg_to_int)
            
            
            ## match if both close enough and same int as the rounded average
            model.addGenConstrAnd(mcn_match[t,s], [mcn_close_to_int[t,s], mcn_close_to_avg[t,s]]) 
            model.addGenConstrAnd(mcn_match_and_avg_at_int[t,s], [mcn_match[t,s], mcn_avg_close_to_int[t,s]]) 

            ## constraint for TCN-based CNA
            model.addConstr((mcn_gain[t,s]==1) >> (mcn_int[t,s] >= mcn_wt_copies + 1))
            model.addConstr((mcn_loss[t,s]==1) >> (mcn_int[t,s] <= mcn_wt_copies - 1))
            model.addGenConstrOr(mcn_cna[t,s], [mcn_gain[t,s], mcn_loss[t,s]])



            # =============================================================================
            # combined TCN and MCN
            # =============================================================================
       
            ## check for both TCN and MCN match
            model.addGenConstrAnd(match_both[t,s], [tcn_match_and_avg_at_int[t,s], mcn_match_and_avg_at_int[t,s]]) 
                
            ## check for CNA in TCN or MCN
            model.addGenConstrOr(is_cna[t,s], [tcn_cna[t,s], mcn_cna[t,s]]) 
                
            ## check if it has homdel
            model.addGenConstrIndicator(is_homdel[t,s], 1, tcn[t,s], GRB.LESS_EQUAL, 0.5)
            
            ## check if segment matches and is large and has a CNA 
            model.addGenConstrAnd(match_both_and_large_enough[t,s], [match_both[t,s], large_enough[t,s]])
            model.addGenConstrAnd(match_both_and_large_enough_and_cna[t,s], [match_both_and_large_enough[t,s], is_cna[t,s]])

    for s in Segments:    
        model.addGenConstrIndicator(allmatch[s], 1, gb.quicksum(match_both_and_large_enough_and_cna[(t, s)] for t in Samples), GRB.GREATER_EQUAL, rho*n_Samples)

    # get total homdel Mb and number of segments with CNAs for each sample
    for t in Samples:
        model.addConstr(homdel_mb[t] == gb.quicksum(dat.loc[(t,s)].mb * is_homdel[(t, s)] for s in Segments))
        model.addConstr(n_cna_segments_per_sample[t] == gb.quicksum(is_cna[(t,s)] for s in Segments))




    # =============================================================================
    # define objectives
    # =============================================================================

    # objective 1: number of segments with clonal SCNAs (the same CNA in 1+ allele, present in rho+ % of samples)
    model.addConstr(n_clonal == gb.quicksum(allmatch[s] for s in Segments))
    
    if(obj2_clonalonly==0):
        # objective 2a, 2b: minimize the combined error among all segments
        model.addConstr(tcn_error_clonal == gb.quicksum(tcn_int_err[t, s] for t in Samples for s in Segments))
        model.addConstr(mcn_error_clonal == gb.quicksum(mcn_int_err[t, s] for t in Samples for s in Segments))    
    
    else:
        # objective 2a, 2b: minimize the combined error among CLONAL segments 
        tcn_int_err_term = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name="tcn_int_err_term")
        mcn_int_err_term = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name="mcn_int_err_term")
        for t in Samples:
            for s in Segments:
                model.addConstr(tcn_int_err_term[t, s] <= tcn_int_err[t, s])
                model.addConstr(tcn_int_err_term[t, s] <= allmatch[s])
                model.addConstr(tcn_int_err_term[t, s] >= tcn_int_err[t, s] - (1 - allmatch[s]))
                model.addConstr(tcn_int_err_term[t, s] >= 0)    
                model.addConstr(mcn_int_err_term[t, s] <= mcn_int_err[t, s])
                model.addConstr(mcn_int_err_term[t, s] <= allmatch[s])
                model.addConstr(mcn_int_err_term[t, s] >= mcn_int_err[t, s] - (1 - allmatch[s]))
                model.addConstr(mcn_int_err_term[t, s] >= 0)              
        model.addConstr(tcn_error_clonal == gb.quicksum(tcn_int_err_term[t, s] for t in Samples for s in Segments))
        model.addConstr(mcn_error_clonal == gb.quicksum(mcn_int_err_term[t, s] for t in Samples for s in Segments))

    # Optimize with stagnation callback  
    model.setObjectiveN(n_clonal, index=0, priority=2, weight=1, name='N clonal segs')
    model.setObjectiveN(-tcn_error_clonal, index=1, priority=1, weight=1-mcn_weight, name='TCN error')
    model.setObjectiveN(-mcn_error_clonal, index=1, priority=1, weight=mcn_weight, name='MCN error')

    model.ModelSense = GRB.MAXIMIZE
    model.update()    
    model.optimize(callback)
    
    #model.computeIIS()
    #model.write("/Users/alexgorelick/Desktop/infeasible_model.ilp")
    return model




