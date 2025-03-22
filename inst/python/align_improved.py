
import pandas as pd
import gurobipy as gp
import numpy as np
from gurobipy import GRB

## non-allele specific version.
## dat should be a data.frame object from R with columns: "sample", "segment", "logR", "GC"
def do_CNalign(dat, min_ploidy=1.7, max_ploidy=6.0, min_purity=0.1, max_purity=1.0, max_homdel_mb=10, t_diff=0.05, t_err=0.1, rho=0.85, min_cna_segments=5, gurobi_license='', epsilon=1e-4):

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
    num_samples = len(Samples)
    num_segments = len(Segments)
    dat.set_index(['sample','segment'], inplace=True) ## set indices: sample, segment

    ## print out message with input parameters 
    print('\n-------------------------------------')
    print('Running optimal alignment with parameters:')
    print('Gurobi license: '+gurobi_license)  
    print('ploidy range: ['+str(min_ploidy)+'-'+str(max_ploidy)+']')  
    print('purity range: ['+str(min_purity)+'-'+str(max_purity)+']')  
    print('Max homdel length allowed: '+str(max_homdel_mb))  
    print('rho (max fraction of samples with matching segment): '+str(rho))  
    print('t_diff (max diff from integer CN): '+str(t_diff))  
    print('t_err (max diff from integer CN): '+str(t_err))  
    print('# samples in file: '+str(num_samples))
    print('# segments in file: '+str(num_segments))
    print(f'Samples: {Samples}')
    print('-------------------------------------')

    env = gp.Env(params=params)
    model = gp.Model(env=env)

    # sample-level
    pl = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='pl', lb=min_ploidy, ub=max_ploidy)
    z = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='z', lb=1/max_purity, ub=1/min_purity)
    homdel_mb = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='homdel_mb', lb=0, ub=max_homdel_mb) ## get number of homdel segments in each sample
    n_cna_segments = model.addVars(Samples, vtype=GRB.INTEGER, name='n_cna_segments', lb=0, ub=num_segments) ## get number of homdel segments in each sample
    segment_pl = model.addVars(Samples, vtype=GRB.CONTINUOUS, name='segment_pl', lb=min_ploidy, ub=max_ploidy)

    # segment-level
    n_bar = model.addVars(Segments, vtype=GRB.CONTINUOUS, name='n_bar')
    allmatch = model.addVars(Segments, vtype=GRB.BINARY, name='allmatch')

    # sample+segment level
    n = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='n') # total copy number cannot be negative
    n_int = model.addVars(Samples, Segments, vtype=GRB.INTEGER, name='n_int')
    dist_from_int = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, lb=0, name='dist_from_int')
    err = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='err')
    diff = model.addVars(Samples, Segments, vtype=GRB.CONTINUOUS, name='diff')
    match_diff = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_diff')
    match_err = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_err')
    match = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match')
    match_with_cna = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='match_with_cna')
    homdel = model.addVars(Samples, Segments, vtype=GRB.BINARY, name='homdel')
    DEL = model.addVars(Samples, Segments, vtype=GRB.BINARY,name='DEL')
    DUP = model.addVars(Samples, Segments, vtype=GRB.BINARY,name='DUP')
    CNA = model.addVars(Samples, Segments, vtype=GRB.BINARY,name='CNA')

    dist = model.addVars(Samples, Samples, Segments, vtype=GRB.CONTINUOUS, name='dist', lb=0)
    L1 = model.addVars(Samples, Samples, vtype=GRB.CONTINUOUS, name='L1', lb=0)
    total_dist_from_int = model.addVar(vtype=GRB.CONTINUOUS, name='total_dist_from_int', lb=0)
    total_tree_distance = model.addVar(vtype=GRB.CONTINUOUS, name='total_tree_distance', lb=0)

    # calculate TCN
    for t in Samples:
        if t=='Diploid':
            # the first sample should be the normal/diploid
            model.addConstr(z['Diploid']==1)
            model.addConstr(pl['Diploid']==2)

        # the segments should approximate the sample's ploidy
        silence = model.addConstr(segment_pl[t] == gp.quicksum(n[t,s]*dat.loc[t,s].mb for s in Segments) / np.sum(dat.loc[t,:].mb))
        silence = model.addConstr(segment_pl[t] >= pl[t] - 0.25) 
        silence = model.addConstr(segment_pl[t] <= pl[t] + 0.25) 

        for s in Segments:
            ## constants
            logr = dat.loc[t,s].logR
            g = dat.loc[t,s].GC 
            r = 2**logr
            silence = model.addConstr(n[t,s] == z[t]*r*2 + r*pl[t] - 2*r - g*z[t] + g) # n=TCN
            silence = model.addConstr(n_int[t,s] >= n[t,s] - 0.5) # n_int is rounded TCN
            silence = model.addConstr(n_int[t,s] <= n[t,s] + 0.5 - epsilon) 
            silence = model.addConstr(dist_from_int[t,s] >= n_int[t,s] - n[t,s])
            silence = model.addConstr(dist_from_int[t,s] >= n[t,s] - n_int[t,s])
            silence = model.addGenConstrIndicator(homdel[t,s], 1, n[t,s], GRB.LESS_EQUAL, 0.5-epsilon) # homdel regions
            silence = model.addGenConstrIndicator(DEL[(t,s)], 1, n[(t, s)], GRB.LESS_EQUAL, 1.5-epsilon) 
            silence = model.addGenConstrIndicator(DUP[(t,s)], 1, n[(t, s)], GRB.GREATER_EQUAL, 2.5) 
            silence = model.addGenConstrOr(CNA[(t, s)], [DEL[(t, s)], DUP[(t, s)]])

    # total length with homdels per sample
    for t in Samples:
        silence = model.addConstr(n_cna_segments[t] == gp.quicksum(CNA[t,s] for s in Segments))
        silence = model.addConstr(homdel_mb[t] == gp.quicksum(homdel[t,s]*dat.loc[t,s].mb for s in Segments))
        if t!='Diploid':
            silence = model.addConstr(n_cna_segments[t] >= min_cna_segments)

    # L1 distance calculation
    for t1 in Samples:
        for t2 in Samples:
            if t1!=t2:
                silence = model.addConstr(dist[t1,t2,s] >= n_int[t1,s] - n_int[t2,s])
                silence = model.addConstr(dist[t1,t2,s] >= n_int[t2,s] - n_int[t1,s])
                #silence = model.addConstr(dist[t1,t2,s] >= n[t1,s] - n[t2,s])
                #silence = model.addConstr(dist[t1,t2,s] >= n[t2,s] - n[t1,s])
                silence = model.addConstr(L1[t1,t2] == gp.quicksum(dist[t1,t2,s] for s in Segments))
            else:
                silence = model.addConstr(dist[t1,t2,s] == 0)
                silence = model.addConstr(L1[t1,t2] == 0)

    # Objective 1: minimize NJ tree distance sum
    model.addConstr(total_tree_distance == 0.5*gp.quicksum(L1[t1, t2] for t1 in Samples for t2 in Samples))


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # minimum spanning tree
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#    n_nodes = num_samples
#
#    # node order in graph
#    order = model.addVars(n_nodes, vtype=GRB.INTEGER, name="order")
#
#    # Binary edge variables: x[i, j] = 1 if there is an edge from node i to j
#    edges = model.addVars(n_nodes, n_nodes, vtype=GRB.BINARY, name="edges")
#
#    # select root node in order
#    model.addConstr(order[0] == 1)
#
#    # No bidirectional edges
#    for i in range(n_nodes):
#        for j in range(n_nodes):
#            model.addConstr(edges[i, j] + edges[j, i] <= 1) 
#
#    # Root node constraints (only outgoing edges)
#    model.addConstr(gp.quicksum(edges[j, 0] for j in range(1, n_nodes)) == 0)  # No incoming edges to root
#    model.addConstr(gp.quicksum(edges[0, j] for j in range(1, n_nodes)) == 1)  # Root must have 1 outgoing edge
#
#    # All samples are nodes
#    for i in range(1,n_nodes):
#        model.addConstr(gp.quicksum(edges[j, i] for j in range(n_nodes) if i != j) == 1) # all nodes have 1 parent
#
#    # Ensure tree connectivity (number of edges must be n_nodes - 1)
#    model.addConstr(gp.quicksum(edges[i, j] for i in range(n_nodes) for j in range(n_nodes)) == n_nodes - 1) # number of edges for a rooted tree
#
#    # Miller-Tucker-Zemlin (MTZ) constraints
#    for i in range(1, n_nodes):  # Skip root
#        for j in range(n_nodes):
#            if i != j:
#                model.addConstr(order[i] >= order[j] + 1 - n_nodes * (1 - edges[i, j]))
#
#    # Objective 1: minimize total weight of selected edges
#    model.addConstr(total_tree_distance == gp.quicksum(L1[Samples[i], Samples[j]] * edges[i, j] for i in range(n_nodes) for j in range(n_nodes)))


    # Objective 2: minimize overall distance from integer copies
    model.addConstr(total_dist_from_int == gp.quicksum(dist_from_int[t, s] for t in Samples for s in Segments) / (num_samples*num_segments) )

    # prioritize best fit, then min evolution
    #model.setObjectiveN(total_dist_from_int, 0, 1, name='Nearest integer')
    #model.setObjectiveN(total_tree_distance, 1, 0, name='Evolutionary distance')
    #model.ModelSense = GRB.MINIMIZE

    # prioritize best fit, then min evolution
    model.setObjectiveN(total_tree_distance, 0, 1, name='Evolutionary distance')
    model.setObjectiveN(total_dist_from_int, 1, 0, name='Nearest integer')
    model.ModelSense = GRB.MINIMIZE

    model.update()
    model.optimize()
    return model


