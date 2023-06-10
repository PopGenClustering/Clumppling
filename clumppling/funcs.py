import numpy as np
import pandas as pd
import os
import shutil
import json
import sys

from itertools import combinations, product 
from collections import defaultdict, Counter
from scipy.spatial.distance import cdist
import cvxpy as cp

import networkx as nx
from TracyWidom import TracyWidom
import community.community_louvain as community_louvain

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import ConnectionPatch


def load_default_parameters(parameters):
        
    bool_params = ['merge_cls','use_rep','cd_default','plot_modes','plot_modes_withinK','plot_major_modes','plot_all_modes']
    for par in bool_params:
        if par in parameters:
            parameters[par] = True if parameters[par] in ['T','t','True','true','Yes','Y','y',1] else False
        
    return parameters


def load_parameters(params_path):

    with open(params_path) as f:
        parameters = json.load(f)
        
    bool_params = ['merge_cls','use_rep','cd_default','plot_modes','plot_modes_withinK','plot_major_modes','plot_all_modes']
    for par in bool_params:
        if par in parameters:
            parameters[par] = True if parameters[par] in ['T','t','True','true','Yes','Y','y',1] else False
        
    return parameters


def display_parameters(input_path,input_format,output_path,params_path,parameters):
    disp = []
    disp.append("=========== Parameters ===========")
    disp.append("----------- [Required] -----------")
    disp.append("Input path: {}".format(input_path))
    disp.append("Input data format: {}".format(input_format))
    disp.append("Output path: {}".format(output_path))
    disp.append("Parameter file path: {}".format(params_path))
    disp.append("------------ [Methods] -----------")
    disp.append("Using default community detection method: {}".format(parameters['cd_default']))
    disp.append("Community detection parameter: {}".format(parameters['cd_parameter']))
    disp.append("Community detection modularity threshold: {}".format(parameters['cd_modularity_threshold']))
    disp.append("Using a representative replicate as the mode consensus: {}".format(parameters['use_rep']))
    disp.append("Merging all possible pairs of clusters when aligning two replicates with K differing by one: {}".format(parameters['merge_cls']))
    disp.append("----------- [Plotting] -----------")
    disp.append("Providing customized colormap: {}".format("False" if parameters['custom_cmap']=='' else parameters['custom_cmap']))
    disp.append("Plotting aligned modes on top of a multipartitie graph: {}".format(parameters['plot_modes']))
    disp.append("Plotting modes of the same K: {}".format(parameters['plot_modes_withinK']))
    disp.append("Plotting major modes: {}".format(parameters['plot_major_modes']))
    disp.append("Plotting all modes: {}".format(parameters['plot_all_modes']))
    disp.append("----------------------------------")
    disp.append("==================================")
    
    return "\n".join(disp)


def load_inputs(data_path,output_path,input_format):

    files = os.listdir(data_path)
    indivq_ext = False
    if input_format == "structure":
        input_files = [i for i in files if i.endswith('_f')]
    elif input_format == "fastStructure":
        input_files = [i for i in files if i.endswith('.meanQ')]
    elif input_format == "admixture" or input_format == "generalQ":
        input_files = [i for i in files if i.endswith('.Q')]
        if len(input_files)==0: # indivq files
            input_files = [i for i in files if i.endswith('.indivq')]
            indivq_ext = True
    else:
        sys.exit("ERROR: Please specify one of the following for input_format: structure, admixture, fastStructure, and generalQ.")

    R = len(input_files)
    if R==0:  # check if no input data of the right form 
        sys.exit("ERROR: No input files detected. Please double check input_path and input_format.")
    
    recode_path = os.path.join(output_path,"input")
    if os.path.exists(recode_path):
        shutil.rmtree(recode_path)
    if not os.path.exists(recode_path):
        os.mkdir(recode_path)
    
    Q_files = list()
    Q_list = list()
    K_list = list()
    
    if input_format =="structure":
        for r in range(R):
            input_file = input_files[r]
            with open(os.path.join(data_path,input_file)) as file:
                lines = file.readlines()
            
            # get only the membership coefficients part
            res_start = lines.index("Inferred ancestry of individuals:\n")
            lines = lines[(res_start+1):]
            res_end = lines.index("\n")
            lines = lines[:res_end]
            
            # extract membership matrix
            Q = [[float(i) for i in l.split()[5:]] for l in lines[1:]]
            Q = np.array(Q)
            K = Q.shape[1]
            Q_list.append(Q)
            K_list.append(K)

    elif indivq_ext: # input_format =="admixture"
        for r in range(R):
            input_file = input_files[r]
            Q = np.loadtxt(os.path.join(data_path,input_file),dtype=str)[:,5:].astype(float)
            K = Q.shape[1]
            Q_list.append(Q)
            K_list.append(K)

    else: # if input_format =="admixture" or "fastStructure" or "generalQ"
        for r in range(R):
            input_file = input_files[r]
            Q = np.loadtxt(os.path.join(data_path,input_file)).astype(float)
            K = Q.shape[1]
            Q_list.append(Q)
            K_list.append(K)
    
    # reorder by consecutive K values   
    sorted_idx = list(np.argsort(K_list))
    K_list = [K_list[i] for i in sorted_idx]
    Q_list = [Q_list[i] for i in sorted_idx]
    reordered_input_files = [input_files[i] for i in sorted_idx]
    
    # summarize input
    R = len(Q_list)
    N = Q_list[0].shape[0]
    K_range = np.sort(np.unique(K_list))
    K_max = max(K_range)
    K2IDs = {k:[] for k in K_range}

    # save to Q files
    with open(os.path.join(recode_path,'input_files.txt'),'w') as file:
        for r in range(R):
            K = K_list[r]
            K2IDs[K].append(r)
            Q_file = '{}_K{}R{}'.format(r+1,K,len(K2IDs[K])) # (ID, K, ID in K)
            np.savetxt(os.path.join(recode_path,"{}.Q".format(Q_file)), Q_list[r], delimiter=' ')
            Q_files.append(Q_file)
            file.write('{},{}\n'.format(Q_file,reordered_input_files[r]))


    # extract individual info if available
    if input_format =="structure":
        df_ind_info = pd.DataFrame([l.split()[0:4] for l in lines[1:]])
        df_ind_info.set_index([0],inplace=True)
        header = lines[0].split()[:3]
        header[2] = header[2].split(":")[0]
        df_ind_info.columns = header
        df_ind_info.to_csv(os.path.join(recode_path,'ind_info.txt'), sep=' ')
        
    if indivq_ext:
        df_ind_info = pd.DataFrame(np.loadtxt(os.path.join(data_path,input_file),dtype=str)[:,:4])
        df_ind_info.set_index([0],inplace=True)
        df_ind_info.columns = np.array(['Label','(%Miss)','Pop'])
        df_ind_info.to_csv(os.path.join(recode_path,'ind_info.txt'), sep=' ')

    return Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs




def align_ILP(P,Q):
    
    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    N_ind = Q.shape[0] #individuals
    assert(K_Q>=K_P) # Q has more clusters
    
    # compute distance
    D = cdist(Q.T, P.T, metric='sqeuclidean')/N_ind # seuclidean, euclidean  
    d = D.flatten()
    
    # constraints
    # sum of row = 1
    A = np.zeros((K_Q,d.shape[0]))
    for i in range(K_Q):
        A[i,i*K_P:(i+1)*K_P] = np.ones((1,K_P))
    A_ones = np.ones((K_Q))
    
    # sum of column >= 1
    B = np.zeros((K_P,d.shape[0]))
    for i in range(K_P):
        B[i,np.arange(i,d.shape[0],K_P)] = np.ones((1,K_Q))
    B_ones = np.ones((K_P))
 
    # ILP
    w = cp.Variable(d.shape[0], boolean = True)
    constraints = [A @ w == A_ones, B @ w >= B_ones]
    obj = d.T @ w
    problem = cp.Problem(cp.Minimize(obj), constraints)
    problem.solve(solver=cp.GLPK_MI) #
    
    # optimal solution
    opt_obj = problem.value
    opt_w = problem.variables()[0].value.astype(int)
    opt_W = np.reshape(opt_w,D.shape)

    # matching index Q_idx to P_idx
    idxQ2P = np.where(opt_W==1)[1]
    

    return opt_obj, idxQ2P#, idxP2Q


# align P and Q by enumerating all combinations of two clusters when K is differ by 1
def align_ILP_diff1(P,Q):

    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    N_ind = Q.shape[0] #individuals
    assert(K_Q-K_P==1)

    # all combinations of Q into the same size of P
    Q_2combs = list(combinations(np.arange(K_Q),2))

    best_obj = np.inf
    idxQ2P = None 

    for pair in Q_2combs:
        
        # remove matrix
        old_idx = np.arange(K_Q)
        old_idx = np.delete(old_idx, [pair[0],pair[1]])
        
        # update matrix
        new_idx = Q.shape[1]-2
        Q_comb = np.hstack([Q[:,old_idx],np.expand_dims(Q[:,pair[0]]+Q[:,pair[1]],1)])

        opt_obj, idxQ2P_comb = align_ILP(P,Q_comb)

        if opt_obj<best_obj:
            best_obj = opt_obj
            idxQ2P = np.zeros(K_Q)

            idxQ2P[old_idx] = idxQ2P_comb[np.arange(len(old_idx))]
            idxQ2P[pair[0]] = idxQ2P_comb[new_idx]
            idxQ2P[pair[1]] = idxQ2P_comb[new_idx]
            
    idxQ2P = idxQ2P.astype(int)

    return best_obj, idxQ2P


def align_ILP_weighted(P,Q,weight):
    
    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    N_ind = Q.shape[0] #individuals
    assert(K_Q>=K_P)

    # compute distance
    dist_elem = (np.expand_dims(Q, axis=2)-np.expand_dims(P, axis=1))**2
    
    w_dist = dist_elem*np.expand_dims(weight,(1,2))
    D = np.sum(w_dist,axis=0)
    d = D.flatten()
    
    # constraints
    # sum of row = 1
    A = np.zeros((K_Q,d.shape[0]))
    for i in range(K_Q):
        A[i,i*K_P:(i+1)*K_P] = np.ones((1,K_P))
    A_ones = np.ones((K_Q))
    
    # sum of column >= 1
    B = np.zeros((K_P,d.shape[0]))
    for i in range(K_P):
        B[i,np.arange(i,d.shape[0],K_P)] = np.ones((1,K_Q))
    B_ones = np.ones((K_P))
 
    # ILP
    w = cp.Variable(d.shape[0], boolean = True)
    constraints = [A @ w == A_ones, B @ w >= B_ones]
    obj = d.T @ w
    problem = cp.Problem(cp.Minimize(obj), constraints)
    problem.solve(solver=cp.GLPK_MI) #
    
    # optimal solution
    opt_obj = problem.value
    opt_w = problem.variables()[0].value.astype(int)
    opt_W = np.reshape(opt_w,D.shape)

    # matching index Q_idx to P_idx
    idxQ2P = np.where(opt_W==1)[1]
    
    return opt_obj, idxQ2P#, idxP2Q


def alignQ_wrtP(P,Q,idxQ2P,merge=True):
    # K1, K2 = P.shape[1], Q.shape[1]
    
    if merge:        
        aligned_Q = np.zeros_like(P)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        idxQ2P = list(idxQ2P)
        dups = np.unique([i for i in idxQ2P if idxQ2P.count(i)>1])

        
        dups_min = defaultdict(lambda: (float('inf'),None))

        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx not in dups:
                aligned_Q[:,p_idx] = Q[:,q_idx]
            else:
                diff = np.linalg.norm(Q[:,q_idx]-P[:,p_idx])
                if dups_min[p_idx][0] > diff:
                    dups_min[p_idx] = (diff,q_idx) 
        
        extra_cnt = 0
        P_dim2 = P.shape[1]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx in dups:
                if q_idx==dups_min[p_idx][1]:
                    aligned_Q[:,p_idx] = Q[:,q_idx]
                else:
                    aligned_Q[:,P_dim2+extra_cnt] = Q[:,q_idx]
                    extra_cnt += 1
            
    return aligned_Q


def cost_membership(P,Q,N_ind):
    return np.linalg.norm(P-Q, ord="fro")**2/(2*N_ind)


def cost_membership_sep(P,Q,idxQ2P):
    N_ind = P.shape[0]
    cost_col = list()
    for p_idx in range(P.shape[1]):
        cost_col.append(list())
    for q_idx in range(Q.shape[1]):
        p_idx = idxQ2P[q_idx]
        cost_col[p_idx].append([np.sum((P[:,p_idx]-Q[:,q_idx])**2)])
    cost = 0    
    for p_idx in range(P.shape[1]):
        if len(cost_col[p_idx])==1:
            cost += np.mean(cost_col[p_idx])
    cost /= 2*N_ind
    return cost


def cost_membership_unmerged(P,Q_aligned):
    assert(P.shape[1]<=Q_aligned.shape[1])
    return cost_membership(P,Q_aligned[:,:P.shape[1]],P.shape[0])


def reverse_alignment(idxQ2P):
    idxP2Q = np.zeros(len(idxQ2P)).astype(int)
    for i,idx in enumerate(idxQ2P):
        idxP2Q[idx] = i
    return idxP2Q


def align_withinK(output_path,Q_list,Q_files,K_range,K2IDs):

    withinK_path = os.path.join(output_path,"alignment_withinK")
    if os.path.exists(withinK_path):
        shutil.rmtree(withinK_path)
    if not os.path.exists(withinK_path):
        os.mkdir(withinK_path)
    
    alignment_withinK = dict()
    cost_withinK = dict()   
        
    for K in K_range:
        with open(os.path.join(withinK_path,"K{}.txt".format(K)),"w") as f:
            f.write('Replicate1-Replicate2,Cost,Alignment\n')
            ids = K2IDs[K]
            n_ids = len(ids)
                
            for i in range(n_ids-1):
                P = Q_list[ids[i]]
                
                for j in range(i+1,n_ids):
                    Q = Q_list[ids[j]]
                    opt_obj, idxQ2P = align_ILP(P, Q)
                    idxP2Q = reverse_alignment(idxQ2P)
                    
                    # calculate cost
                    aligned_Q = alignQ_wrtP(P,Q,idxQ2P,merge=True)
                    cost = cost_membership(P,aligned_Q,P.shape[0])
                    
                    # write output
                    f.write("{}-{},{},{}\n".format(Q_files[ids[i]],Q_files[ids[j]],cost," ".join([str(id+1) for id in idxQ2P])))
                    f.write("{}-{},{},{}\n".format(Q_files[ids[j]],Q_files[ids[i]],cost," ".join([str(id+1) for id in idxP2Q])))
                    
                    alignment_withinK["{}-{}".format(Q_files[ids[i]],Q_files[ids[j]])] = idxQ2P
                    alignment_withinK["{}-{}".format(Q_files[ids[j]],Q_files[ids[i]])] = idxP2Q
                    cost_withinK["{}-{}".format(Q_files[ids[i]],Q_files[ids[j]])] = cost
                    cost_withinK["{}-{}".format(Q_files[ids[j]],Q_files[ids[i]])] = cost

    return alignment_withinK, cost_withinK



def get_adj_mat(cost_mat):
    adj_mat = 1-(cost_mat-np.nanmin(cost_mat))/(np.nanmax(cost_mat)-np.nanmin(cost_mat))
    # adj_mat = 1-(cost_mat-np.nanmin(cost_mat))/(np.nanmax(cost_mat)-np.nanmin(cost_mat))*0.9
    adj_mat = np.nan_to_num(adj_mat,copy=True,nan=0)
    np.fill_diagonal(adj_mat, 1)
    assert(np.mean(np.diag(adj_mat))==1)
    # adj_mat = adj_mat + np.diag(np.ones(adj_mat.shape[0]))
    return adj_mat


def standardize_matrix(W):
    off_diag_idx = np.where(~np.eye(W.shape[0],dtype=bool))
    off_diag_w = W[off_diag_idx]
    W_standardized = np.zeros(W.shape)
    W_standardized[off_diag_idx] = (off_diag_w-off_diag_w.mean())/off_diag_w.std()
    return W_standardized

def normalize_matrix(W):
    return W/np.sqrt(W.shape[0])

def exponentiate_matrix(W,t):
    off_diag_idx = np.where(~np.eye(W.shape[0],dtype=bool))
    W_exp = np.zeros(W.shape)
    W_exp[off_diag_idx] = np.exp(W[off_diag_idx]*t)
    return W_exp

def test_comm_struc(W, alpha = 0.05):
    """Test for the existence of community structure in the graph (implementing Tokuda 2018)

    Parameters
    ----------
    W : edge-weight matrix of the graph
        The adjacency matrix of the the similarity network of replicates, with edges weighted by similarity after optimal alignment
    alpha : the p-value threshold for the test

    Returns
    -------
    has_comm_struc
        a boolean indicating whether there exists community structure in the graph
        True if the null hypothesis is rejected; False if there is no community structure.
    """
    
    # standardization mapping
    W_standardized = standardize_matrix(W)
    T = normalize_matrix(W_standardized)
    
    t = 1/2
    W_exp = exponentiate_matrix(W,t)
    Te = normalize_matrix(standardize_matrix(W_exp))
    
    # Tracy Widom CI
    x = np.linspace(-10, 10, 1001)
    tw = TracyWidom(beta=1)  # allowed beta values are 1, 2, and 4
    cdf = tw.cdf(x)
    
    CI_max = x[np.where(cdf>(1-alpha/4))[0][0]]
    CI_min = x[np.where(cdf<alpha/4)[0][-1]]
    CI_p_max, CI_p_min = CI_max, CI_min
    
    # test
    eigval, eigvec = np.linalg.eig(T)
    eig_T_max, eig_T_min = eigval.max(), eigval.min()
    eigval, eigvec = np.linalg.eig(Te)
    eig_Te_max, eig_Te_min = eigval.max(), eigval.min()
    
    s = int(eig_T_max<CI_max) + int(eig_T_min>CI_min) + int(eig_Te_max<CI_p_max) + int(eig_Te_min>CI_p_min)
    has_comm_struc = s!=4
    
    return has_comm_struc


def reorder_partition_map(partition_map,G):
    pm_keys = np.array(list(partition_map.keys()))
    pm_vals = np.array(list(partition_map.values()))

    modes_by_size = Counter(pm_vals).most_common()
    to_sort = list()
    for i,mbs in enumerate(modes_by_size):
        mode = mbs[0]
        node_indices = pm_keys[np.where(pm_vals==mode)[0]]
        # break tie based on within-mode similarity
        subG = G.subgraph(node_indices)
        subG_weight = list(nx.get_edge_attributes(subG,'weight').values())
        subG_avg_weight = np.mean(subG_weight)
        to_sort.append((mode,mbs[1],subG_avg_weight))
    sorted_list = sorted(to_sort, key=lambda x: (x[1], x[2]), reverse=True)
    reindex_modes = {mode[0]:i for i,mode in enumerate(sorted_list)}
    
    for k in partition_map:
        partition_map[k] = reindex_modes[partition_map[k]]
    return partition_map

def cd_default(G,res=1.05):
    """Community detection using the Louvain method

    Parameters
    ----------
    G : networkx.Graph
        the similarity network of replicates, with edges weighted by similarity after optimal alignment

    Returns
    -------
    partition_map
        a dictionary where keys are the indices of replicates and values are the indices of the communities they belong to
    """

    resolution = res 
    partition_map = community_louvain.best_partition(G,resolution=resolution,random_state=6)
    return partition_map

def cd_custom(G):
    """Customized community detection method (need to be modified)

    Parameters
    ----------
    G : networkx.Graph
        the similarity network of replicates, with edges weighted by similarity after optimal alignment

    Returns
    -------
    partition_map
        a dictionary where keys are the indices of replicates and values are the indices of the communities they belong to
    
    Example: 
        # using mcode for cummunity detection
        # need to install the cdlib package (https://cdlib.readthedocs.io/en/latest/overview.html) and do "from cdlib import algorithms" 
        def cd_custom(G):
            coms = algorithms.mcode(G,"weight")
            partition = coms.communities
            partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
            return partition_map
    """

    # Please comment out the following line and customize your community detection method here.
    # partition_map = {i:0 for i in range(G.number_of_nodes())} 
    from cdlib import algorithms
    coms = algorithms.markov_clustering(G,pruning_threshold=0.2) # scipy 1.8.0
    partition = coms.communities
    partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
    
    return partition_map


def detect_modes(cost_withinK,Q_files,K_range,K2IDs,default_cd,cd_mod_thre,cd_param=1.05):

    modes_allK = dict()
    cost_matrices = dict()
    msg = []
            
    for K in K_range:
    
        # normalize cost matrix and obtain adjacency matrix
        ids = K2IDs[K]
        cost_mat = np.zeros((len(ids),len(ids)))

        for ii,i in enumerate(ids):
            for jj,j in enumerate(ids):
                if i!=j:
                    cost_mat[ii,jj] = cost_withinK["{}-{}".format(Q_files[i],Q_files[j])]
        cost_matrices[K] = cost_mat
        
        n_nodes = cost_mat.shape[0]
        if np.nanmean(cost_mat)<1e-6:
            # adj_mat = np.zeros(cost_mat.shape)
            partition_map = {i:0 for i in range(n_nodes)} 
            msg.append("K={}: average pairwise cost close to 0 -> set to single mode".format(K))          
        else:
            adj_mat = get_adj_mat(cost_mat)

            # test for community structure
            has_comm_struc = test_comm_struc(adj_mat, alpha = 0.01)
            if not has_comm_struc:
                partition_map = {i:0 for i in range(n_nodes)} 
                msg.append("K={}: no significant community structure -> set to single mode".format(K))

            else: # having community structure
                # create graph
                G = nx.from_numpy_matrix(adj_mat)   
                G.remove_edges_from(nx.selfloop_edges(G))
                pos = nx.spring_layout(G)
                
                ########################################
                # community detection to find modes
                ########################################
                if np.nanmean(cost_mat)<1e-6:
                    partition_map = {i:0 for i in range(n_nodes)} 
                else:
                    if default_cd:
                        partition_map = cd_default(G,cd_param)
                    else:
                        partition_map = cd_custom(G)
                    
                    if cd_mod_thre!=-1:
                        mod_res = community_louvain.modularity(partition_map, G)
                        if mod_res<cd_mod_thre: # quality of community detection is low --> no community structure
                            msg.append("K={}: low community detection quality (below modularity threshold) -> set to single mode".format(K))
                            partition_map = {i:0 for i in range(G.number_of_nodes())} 

                    partition_map = reorder_partition_map(partition_map, G)
        ########################################
        # separate modes
        ########################################
        
        modes = defaultdict(list)
        for r in partition_map.keys():
            modes[partition_map[r]].append(r)
        
        modes_allK[K] = modes

    if len(msg)>0:
        msg = "\n".join(msg)
    else:
        msg = ""
            
    return modes_allK, cost_matrices, msg



def C2Gprime(C):
    return 1-np.sqrt(C)

def avg_tril(A):
    return np.sum(np.tril(A, k=-1))/(A.shape[0]*(A.shape[0]-1)/2)


def extract_modes(Q_list,Q_files,modes_allK,alignment_withinK,cost_matrices,K_range,K2IDs, output_path):
 
    rep_modes = defaultdict(list)
    repQ_modes = defaultdict(list)
    avgQ_modes = defaultdict(list)
    alignment_to_modes = defaultdict(list)
    stats = defaultdict(list)
    
    ########################################
    # get results of alignment for each K
    ########################################
    
    for K in K_range:
    
        modes = modes_allK[K]
        cost_mat = cost_matrices[K]


        for mode_idx in range(len(modes.keys())):

            ########################################
            # find the concensus membership in mode
            ########################################
            all_indices = modes[mode_idx]
            if len(all_indices)==1:
                idx = all_indices[0]
                rep_modes[K].append(idx)
                stats[K].append({"size":1,"cost":0,"perf":1})
                minC_rep_idx = idx+K2IDs[K][0]
                P = Q_list[minC_rep_idx]
                repQ_modes[K].append(P)
                avgQ_modes[K].append(P)
                
                mode_alignment= {Q_files[minC_rep_idx]: np.arange(P.shape[1])}
                alignment_to_modes[K].append(mode_alignment)
                
            else:
                comm_cost_mat = cost_mat[np.array(all_indices),:][:,np.array(all_indices)]
                
                # find the one with max row sum (min cost)
                minC_idx = np.argmin(comm_cost_mat.sum(axis=0))
                # minC_self_idx = all_indices[minC_idx]
                minC_rep_idx = all_indices[minC_idx]+K2IDs[K][0]
                rep_modes[K].append(minC_rep_idx)
                # mds_indices.append(all_indices[minC_idx])

                # other non-representative ones
                # other_self_indices = [m for im,m in enumerate(all_indices) if im!=minC_idx]
                other_rep_indices = [m+K2IDs[K][0] for im,m in enumerate(all_indices) if im!=minC_idx]
                
                ########################################
                # compute the average cost/similarity
                ########################################
                
                comm_cost_mat = cost_mat[np.array(all_indices),:][:,np.array(all_indices)]
                Gprime_mat = C2Gprime(comm_cost_mat)
                stats[K].append({"size":len(all_indices),"cost":avg_tril(comm_cost_mat),"perf":avg_tril(Gprime_mat)})
            
                ########################################
                # get average Q over all in the mode
                ########################################
                mode_alignment = dict()
                
                P = Q_list[minC_rep_idx]
                repQ_modes[K].append(P)
                mode_alignment[Q_files[minC_rep_idx]] = np.arange(P.shape[1])
                
                Q_sum = np.zeros_like(P)
                Q_sum += Q_list[minC_rep_idx]
                # Q_aligned_list = list()
                for i,r in enumerate(other_rep_indices):
                    Q = Q_list[r]
                    aligned_idxQ2P = alignment_withinK["{}-{}".format(Q_files[minC_rep_idx],Q_files[r])]
                    mode_alignment[Q_files[r]] = aligned_idxQ2P
                    aligned_Q = alignQ_wrtP(P,Q,aligned_idxQ2P,merge=True)

                    Q_sum += aligned_Q

                avgQ_modes[K].append(Q_sum/len(all_indices))
                
                alignment_to_modes[K].append(mode_alignment)

    # write output
    modes_path = os.path.join(output_path,"modes")
    if os.path.exists(modes_path):
        shutil.rmtree(modes_path)
    if not os.path.exists(modes_path):
        os.mkdir(modes_path)
    
    f1 = open(os.path.join(modes_path,"mode_stats.txt"),"w")    
    f1.write('Mode,Representative,Size,Cost,Performance\n')
    f2 = open(os.path.join(modes_path,"mode_alignments.txt"),"w")  
    f2.write('Mode,Representative,Replicate,Alignment\n')
    
    mode_labels = defaultdict(list)
    
    for K in K_range:
        for mode_idx in range(len(modes_allK[K].keys())):
            mode_label = "K{}M{}".format(K,mode_idx+1)
            mode_labels[K].append(mode_label)
            np.savetxt(os.path.join(modes_path,"{}_rep.Q".format(mode_label)), repQ_modes[K][mode_idx], delimiter=' ')
            np.savetxt(os.path.join(modes_path,"{}_avg.Q".format(mode_label)), avgQ_modes[K][mode_idx], delimiter=' ')
            
            rep = Q_files[rep_modes[K][mode_idx]]
            f1.write('{},{},{},{},{}\n'.format(mode_label,rep,stats[K][mode_idx]["size"],stats[K][mode_idx]["cost"],stats[K][mode_idx]["perf"]))
            for replicate in alignment_to_modes[K][mode_idx]:
                f2.write('{},{},{},{}\n'.format(mode_label,rep,replicate," ".join([str(id+1) for id in alignment_to_modes[K][mode_idx][replicate]])))

    f1.close()
    f2.close()

    # average statistics
    with open(os.path.join(modes_path,"mode_average_stats.txt"),"w") as f:
        f.write('K,Size,NS-Size,Cost,Performance\n')
        for K in K_range:
            total_cost = 0
            total_perf = 0
            non_singleton_s = 0
            for mode_idx in range(len(modes_allK[K].keys())):
                s = stats[K][mode_idx]["size"]
                if s>1:
                    total_cost += stats[K][mode_idx]['cost']*s
                    total_perf += stats[K][mode_idx]['perf']*s
                    non_singleton_s += s
            # avg_cost += s/len(K2IDs[K])*stats[K][mode_idx]['cost']
            # avg_perf += s/len(K2IDs[K])*stats[K][mode_idx]['perf']
            f.write('{},{},{},{},{}\n'.format(K,len(K2IDs[K]),non_singleton_s,total_cost/non_singleton_s,total_perf/non_singleton_s))

    return mode_labels,rep_modes,repQ_modes,avgQ_modes,alignment_to_modes,stats  


def align_ILP_modes_acrossK(consensusQ_modes,mode_labels,K_range,acrossK_path,cons_suffix,merge=False):
    
    # if os.path.exists(acrossK_path):
    #     shutil.rmtree(acrossK_path)
    if not os.path.exists(acrossK_path):
        os.mkdir(acrossK_path)
    
    K_range_sorted = sorted(K_range,reverse=True)
    K_comb = list([(K_range_sorted[i],K_range_sorted[i+1]) for i in range(len(K_range_sorted)-1)])
    K_comb.extend([(K_range_sorted[i],K_range_sorted[i]) for i in range(len(K_range_sorted))])
    
    alignment_acrossK = dict()
    cost_acrossK = dict()
    
    f = open(os.path.join(acrossK_path,"alignment_acrossK_{}.txt".format(cons_suffix)),"w")
    f.write("Mode1-Mode2,Cost,Alignment\n")
    
    best_alignments = dict()
    
    for K1,K2 in K_comb:
        
        best_alignment_idx = None
        best_alignment_obj = np.inf

        rijs = list(product(range(len(consensusQ_modes[K1])),range(len(consensusQ_modes[K2]))))
        
        for i,(ri,rj) in enumerate(rijs):
            
            Q = consensusQ_modes[K1][ri] # Q is the one with more clusters
            P = consensusQ_modes[K2][rj]

            if (K1==K2 and rj==ri):
                # identical alignment (dummy)
                opt_obj = 0
                idxQ2P = np.arange(K1)
            elif merge and (K1-K2)==1:
                opt_obj, idxQ2P = align_ILP_diff1(P, Q)
            else:
                opt_obj, idxQ2P = align_ILP(P, Q)

            pair_label = "{}-{}".format(mode_labels[K2][rj],mode_labels[K1][ri])
            if best_alignment_obj>opt_obj:
                best_alignment_obj = opt_obj
                best_alignment_idx = i
            
            # calculate cost
            aligned_Q = alignQ_wrtP(P,Q,idxQ2P,merge=True)
            cost = cost_membership(P,aligned_Q,P.shape[0])
            f.write("{}-{},{},{}\n".format(mode_labels[K2][rj],mode_labels[K1][ri],cost," ".join([str(id+1) for id in idxQ2P])))

            cost_acrossK[pair_label] = cost
            alignment_acrossK[pair_label] = idxQ2P
            
        best_alignments[(K1,K2)] = rijs[best_alignment_idx]

    f.close()
    
    # best alignments
    best_acrossK = list()
    f = open(os.path.join(acrossK_path,"best_pairs_acrossK_{}.txt".format(cons_suffix)),"w")
    stats = ['cost_merge', 'perf_merge','cost_sep', 'perf_sep']
    f.write('Best Pair, Cost, Performance, Separate-Cluster Cost, Separate-Cluster Performance\n')
    
    for i_K1 in range(len(K_range_sorted)-1):
        K1 = K_range_sorted[i_K1]
        K2 = K_range_sorted[i_K1+1]
        ri,rj = best_alignments[(K1,K2)]
        bali = "{}-{}".format(mode_labels[K2][rj],mode_labels[K1][ri])
        Q = consensusQ_modes[K1][ri]
        P = consensusQ_modes[K2][rj]
        cost_bali = cost_acrossK[bali]
        best_acrossK.append(bali)

        # Q_aligned_unmerged = alignQ_wrtP(P,Q,alignment_acrossK[bali],merge=False)
        # cost_sep = cost_membership(P,aligned_Q,P.shape[0]) #cost_membership(P,Q_aligned_unmerged[:,:P.shape[1]],P.shape[0])
        cost_sep = cost_membership_sep(P,Q,alignment_acrossK[bali])
        
        f.write("{},{},{},{},{}\n".format(bali,cost_bali,C2Gprime(cost_bali),cost_sep,C2Gprime(cost_sep)))
    
    f.close()
    
    return alignment_acrossK, cost_acrossK, best_acrossK


def plot_structure_on_multipartite(K_range,mode_labels,stats,consensusQ_modes,alignment_acrossK,cost_acrossK_cons,best_acrossK,cons_suffix,output_path,plot_flag,cmap):
    mode_numbers = [len(mode_labels[K]) for K in K_range]

    n_row = len(K_range)
    n_col = max(mode_numbers)

    # use the best alignment mode as the anchors to align all the modes
    best_acrossK = best_acrossK[::-1]
    mode2fig_idx = dict()
    for i_K,K in enumerate(K_range): 
        for i_lb,lb in enumerate(mode_labels[K]):
            if n_col==1:
                mode2fig_idx[lb] = i_K
            else:
                mode2fig_idx[lb] = (i_K,i_lb)

    modeQ_path = os.path.join(output_path,"modes_aligned")
    fig_path = os.path.join(output_path,"visualization")
    if not os.path.exists(modeQ_path):
        os.mkdir(modeQ_path)
        
    all_modes_alignment = {lb:[] for K in K_range for lb in mode_labels[K]}
    base_patterns = dict()
    
    if plot_flag:
        fig, axes = plt.subplots(n_row,n_col,figsize=(5*n_col,2*n_row),facecolor='white')
    
    K = K_range[0]
    m1 = best_acrossK[0].split("-")[0]
    assert(m1.split("M")[0].strip("K")==str(K))
    m_m1 = int(m1.split("M")[1])-1
    K_max = max(K_range)
    base_Q = consensusQ_modes[K][m_m1]
    base_patterns[m1] = [i for i in range(K)]
    if plot_flag:
        ax = axes[mode2fig_idx[m1]]
        plot_membership(ax,base_Q,K_max,cmap,"")   
    all_modes_alignment[m1] = [i+1 for i in range(base_Q.shape[1])]
    np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m1,cons_suffix)), base_Q, delimiter=' ')
    
    modes = range(len(consensusQ_modes[K]))
    if len(modes)>1:
        for m in modes:
            if m!=m_m1:
                m2 = "K{}M{}".format(K,m+1)
                # retrieve alignment pattern
                ali_pat = alignment_acrossK["{}-{}".format(m1,m2)]
                Q = consensusQ_modes[K][m]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                base_patterns[m2] = ali_pat
                if plot_flag:
                    ax = axes[mode2fig_idx[m2]]
                    plot_membership(ax,aligned_Q,K_max,cmap,"")
                all_modes_alignment[m2] = [i+1 for i in ali_pat]
                np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m2,cons_suffix)), aligned_Q, delimiter=' ')
    
    for best_pair in best_acrossK:
        m1 = best_pair.split("-")[0]
        m2 = best_pair.split("-")[1]
        
        K = int(m2.split("M")[0].strip("K"))
        
        modes = range(len(consensusQ_modes[K]))
        m_m2 = int(m2.split("M")[1])-1
        Q = consensusQ_modes[K][m_m2]
        m1_K = int(m1.split("M")[0].strip("K"))
        P = consensusQ_modes[m1_K][int(m1.split("M")[1])-1]
        pattern = alignment_acrossK[best_pair]
        aligned_Q, new_pattern = alignQ_wrtP_pattern(P,Q,pattern,merge=False)

        pat = [base_patterns[m1][i] if i<m1_K else i for i in new_pattern]
        aligned_Q = np.zeros_like(aligned_Q)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
        base_patterns[m2] = pat
        
        if plot_flag:
            ax = axes[mode2fig_idx[m2]]     
            plot_membership(ax,aligned_Q,K_max,cmap,"")
        all_modes_alignment[m2] = [i+1 for i in pat]
        np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m2,cons_suffix)), aligned_Q, delimiter=' ')
        
        if len(modes)>1:
            for m in modes:
                if m!=m_m2:
                    m3 = "K{}M{}".format(K,m+1)
                    ali_pat = alignment_acrossK["{}-{}".format(m2,m3)]
                    Q = consensusQ_modes[K][m]
                    P = consensusQ_modes[K][m_m2]
                    aligned_Q, new_pattern = alignQ_wrtP_pattern(P,Q,ali_pat,merge=False)
                    
                    pat = [base_patterns[m2][i] for i in new_pattern ]
                    aligned_Q = np.zeros_like(aligned_Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
                    base_patterns[m3] = pat
                    
                    if plot_flag:
                        ax = axes[mode2fig_idx[m3]]  
                        plot_membership(ax,aligned_Q,K_max,cmap,"")
                    all_modes_alignment[m3] = [i+1 for i in pat]
                    np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m3,cons_suffix)), aligned_Q, delimiter=' ')            
    
    # write the alignment pattern to result file
    with open(os.path.join(modeQ_path,"all_modes_alignment_{}.txt".format(cons_suffix)),'w') as f:
        for k in all_modes_alignment:
            f.write('{}:{}\n'.format(k," ".join([str(i) for i in all_modes_alignment[k]])))


    if plot_flag:
        # polish the figure
        if n_col>1:
            for j in range(n_col):
                axes[0,j].set_title("Mode {}".format(j+1), fontsize=18, y=1.05)
            for i in range(n_row):
                axes[i,0].set_ylabel("K={}".format(K_range[i]), rotation=0, fontsize=18, labelpad=30, va="center")
        else:
            axes[0].set_title("Mode 1", fontsize=18, y=1.05)
            for i in range(n_row):
                axes[i].set_ylabel("K={}".format(K_range[i]), rotation=0, fontsize=18, labelpad=30, va="center")
        
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=1.2)
        textbox_props = dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='none', pad=0.1)
        
        def get_center_coord(coord1,coord2,frac):
          return (coord2-coord1)*frac+coord1

        if n_col>1:
            for i in range(n_row):
                for j in range(n_col):
                    ax = axes[i,j]
                    ax.set_xticks([])
                    if j>0:
                        ax.set_yticks([])
                    else:
                        ax.set_yticks([0,0.5,1])     
                    if (i,j) not in mode2fig_idx.values():
                        for pos in ['right', 'top', 'bottom', 'left']:
                            ax.spines[pos].set_visible(False)  

            # add edges
            transFigure = fig.transFigure.inverted()
            cost_max = max(cost_acrossK_cons.values())
            cost_min = min(cost_acrossK_cons.values())
            N = Q.shape[0]
            for i in range(n_row-1):
                for j in range(mode_numbers[i]):
                    for jj in range(mode_numbers[i+1]):

                        k1k2_pair = "{}-{}".format(mode_labels[K_range[i]][j],mode_labels[K_range[i+1]][jj])
                        cost = float(cost_acrossK_cons[k1k2_pair])
                        edge_w = 0.9-(cost-cost_min)/(cost_max-cost_min)*0.7
                        
                        cpatch = ConnectionPatch((.5, 0), (.5, 1), "axes fraction", "axes fraction",
                          axesA=axes[i,j], axesB=axes[i+1,jj], zorder=0, 
                          color=plt.cm.Greys(edge_w), linewidth=5, alpha=0.8)
                        axes[0,0].add_artist(cpatch)

                        coord1 = transFigure.transform(axes[i,j].transData.transform([N//2,-0.05]))
                        coord2 = transFigure.transform(axes[i+1,jj].transData.transform([N//2,1.05]))

                        text_frac = 0.25+(jj%2)*0.5 # make the text sparse when edges are dense
                        plt.text(get_center_coord(coord1[0],coord2[0],text_frac),get_center_coord(coord1[1],coord2[1],text_frac), 
                            "{:.3f}".format(cost), fontsize=14, transform=fig.transFigure, bbox=textbox_props, 
                            zorder=9, verticalalignment='center',horizontalalignment='center')
        else:
            transFigure = fig.transFigure.inverted()
            cost_max = max(cost_acrossK_cons.values())
            cost_min = min(cost_acrossK_cons.values())
            N = Q.shape[0]

            for i in range(n_row-1):
                k1k2_pair = "{}-{}".format(mode_labels[K_range[i]][0],mode_labels[K_range[i+1]][0])
                cost = float(cost_acrossK_cons[k1k2_pair])
                edge_w = 0.9-(cost-cost_min)/(cost_max-cost_min)*0.7
                
                cpatch = ConnectionPatch((.5, 0), (.5, 1), "axes fraction", "axes fraction",
                  axesA=axes[i], axesB=axes[i+1], zorder=0, 
                  color=plt.cm.Greys(edge_w), linewidth=5, alpha=0.8)
                axes[0].add_artist(cpatch)
                coord1 = transFigure.transform(axes[i].transData.transform([N//2,-0.05]))
                coord2 = transFigure.transform(axes[i+1].transData.transform([N//2,1.05]))
                text_frac = 0.5 # make the text sparse when edges are dense
                plt.text(get_center_coord(coord1[0],coord2[0],text_frac),get_center_coord(coord1[1],coord2[1],text_frac), 
                    "{:.3f}".format(cost), fontsize=14, transform=fig.transFigure, bbox=textbox_props, 
                    zorder=9, verticalalignment='center',horizontalalignment='center')


        fig.savefig(os.path.join(fig_path,"modes_aligned_multipartite_{}.png".format(cons_suffix)), bbox_inches='tight',dpi=300)
        plt.close(fig)
    

# Just for CV data
def plot_structure_on_multipartite_manuscript(K_range,mode_labels,stats,consensusQ_modes,alignment_acrossK,cost_acrossK_cons,best_acrossK,cons_suffix,output_path,plot_flag,cmap):
    

    mode_numbers = [len(mode_labels[K]) for K in K_range]
    if plot_flag:
        n_row = len(K_range)
        n_col = max(mode_numbers)

    # use the best alignment mode as the anchors to align all the modes
    best_acrossK = best_acrossK[::-1]
    mode2fig_idx = dict()
    for i_K,K in enumerate(K_range): 
        for i_lb,lb in enumerate(mode_labels[K]):
            if n_col==1:
                mode2fig_idx[lb] = i_K
            else:
                mode2fig_idx[lb] = (i_K,i_lb)

    modeQ_path = os.path.join(output_path,"modes_aligned")
    fig_path = os.path.join(output_path,"visualization")
    if not os.path.exists(modeQ_path):
        os.mkdir(modeQ_path)
        
    all_modes_alignment = {lb:[] for K in K_range for lb in mode_labels[K]}
    base_patterns = dict()
    
    if plot_flag:
        fig, axes = plt.subplots(n_row,n_col,figsize=(5*n_col,2*n_row),facecolor='white')
    
    K = K_range[0]
    m1 = best_acrossK[0].split("-")[0]
    assert(m1.split("M")[0].strip("K")==str(K))
    m_m1 = int(m1.split("M")[1])-1
    K_max = max(K_range)
    base_Q = consensusQ_modes[K][m_m1]
    base_patterns[m1] = [i for i in range(K)]
    if plot_flag:
        ax = axes[mode2fig_idx[m1]]
        plot_membership(ax,base_Q,K_max,cmap,"") 

    all_modes_alignment[m1] = [i+1 for i in range(base_Q.shape[1])]
    np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m1,cons_suffix)), base_Q, delimiter=' ')
    
    modes = range(len(consensusQ_modes[K]))
    if len(modes)>1:
        for m in modes:
            if m!=m_m1:
                m2 = "K{}M{}".format(K,m+1)
                # retrieve alignment pattern
                ali_pat = alignment_acrossK["{}-{}".format(m1,m2)]
                Q = consensusQ_modes[K][m]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                base_patterns[m2] = ali_pat
                if plot_flag:
                    ax = axes[mode2fig_idx[m2]]
                    # plot_membership(ax,aligned_Q,K_max,cmap,"")
                    if m2 in ["K5M1","K5M2","K5M3"]:    
                        plot_membership(ax,aligned_Q[:,[3,0,4,1,2]],K_max,cm.colors.ListedColormap(cmap([3,0,4,1,2])),"") 
                    else:
                        plot_membership(ax,aligned_Q,K_max,cmap,"") 
                all_modes_alignment[m2] = [i+1 for i in ali_pat]
                np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m2,cons_suffix)), aligned_Q, delimiter=' ')
    
    for best_pair in best_acrossK:
        m1 = best_pair.split("-")[0]
        m2 = best_pair.split("-")[1]
        
        K = int(m2.split("M")[0].strip("K"))
        
        modes = range(len(consensusQ_modes[K]))
        m_m2 = int(m2.split("M")[1])-1
        Q = consensusQ_modes[K][m_m2]
        m1_K = int(m1.split("M")[0].strip("K"))
        P = consensusQ_modes[m1_K][int(m1.split("M")[1])-1]
        pattern = alignment_acrossK[best_pair]
        aligned_Q, new_pattern = alignQ_wrtP_pattern(P,Q,pattern,merge=False)

        pat = [base_patterns[m1][i] if i<m1_K else i for i in new_pattern]
        aligned_Q = np.zeros_like(aligned_Q)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
        base_patterns[m2] = pat
        
        if plot_flag:
            ax = axes[mode2fig_idx[m2]]     
            # plot_membership(ax,aligned_Q,K_max,cmap,"")
            if m2 in ["K5M1","K5M2","K5M3"]:    
                plot_membership(ax,aligned_Q[:,[3,0,4,1,2]],K_max,cm.colors.ListedColormap(cmap([3,0,4,1,2])),"") 
            else:
                plot_membership(ax,aligned_Q,K_max,cmap,"") 
        all_modes_alignment[m2] = [i+1 for i in pat]
        np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m2,cons_suffix)), aligned_Q, delimiter=' ')
        
        if len(modes)>1:
            for m in modes:
                if m!=m_m2:
                    m3 = "K{}M{}".format(K,m+1)
                    ali_pat = alignment_acrossK["{}-{}".format(m2,m3)]
                    Q = consensusQ_modes[K][m]
                    P = consensusQ_modes[K][m_m2]
                    aligned_Q, new_pattern = alignQ_wrtP_pattern(P,Q,ali_pat,merge=False)
                    
                    pat = [base_patterns[m2][i] for i in new_pattern ]
                    aligned_Q = np.zeros_like(aligned_Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
                    base_patterns[m3] = pat
                    
                    if plot_flag:
                        ax = axes[mode2fig_idx[m3]]  
                        # plot_membership(ax,aligned_Q,K_max,cmap,"")
                        if m3=="K4M1":
                            plot_membership(ax,aligned_Q[:,[3,0,1,2]],K_max,cm.colors.ListedColormap(cmap([3,0,1,2,4])),"")
                        elif m3 in ["K5M1","K5M2","K5M3"]:    
                            plot_membership(ax,aligned_Q[:,[3,0,4,1,2]],K_max,cm.colors.ListedColormap(cmap([3,0,4,1,2])),"") 
                        else:
                            plot_membership(ax,aligned_Q,K_max,cmap,"") 
                    all_modes_alignment[m3] = [i+1 for i in pat]
                    np.savetxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m3,cons_suffix)), aligned_Q, delimiter=' ')            
    
    # write the alignment pattern to result file
    with open(os.path.join(modeQ_path,"all_modes_alignment_{}.txt".format(cons_suffix)),'w') as f:
        for k in all_modes_alignment:
            f.write('{}:{}\n'.format(k," ".join([str(i) for i in all_modes_alignment[k]])))


    if plot_flag:
        # polish the figure
        if n_col>1:
            for j in range(n_col):
                axes[0,j].set_title("Mode {}".format(j+1), fontsize=22, y=1.05)
            for i in range(n_row):
                axes[i,0].set_ylabel("K={}".format(K_range[i]), rotation=0, fontsize=22, labelpad=30, va="center")
        else:
            axes[0].set_title("Mode 1", fontsize=22, y=1.05)
            for i in range(n_row):
                axes[i].set_ylabel("K={}".format(K_range[i]), rotation=0, fontsize=22, labelpad=30, va="center")
        
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=2)
        textbox_props = dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='none', pad=0.1)
        
        def get_center_coord(coord1,coord2,frac):
          return (coord2-coord1)*frac+coord1

        if n_col>1:
            for i in range(n_row):
                for j in range(n_col):
                    ax = axes[i,j]
                    ax.set_xticks([])
                    if j>0:
                        ax.set_yticks([])
                    else:
                        ax.set_yticks([0,0.5,1]) 
                        ax.tick_params(axis='y', labelsize=18)    
                    if (i,j) not in mode2fig_idx.values():
                        for pos in ['right', 'top', 'bottom', 'left']:
                            ax.spines[pos].set_visible(False)  

            # add edges
            transFigure = fig.transFigure.inverted()
            cost_max = max(cost_acrossK_cons.values())
            cost_min = min(cost_acrossK_cons.values())
            N = Q.shape[0]
            for i in range(n_row-1):
                for j in range(mode_numbers[i]):
                    for jj in range(mode_numbers[i+1]):

                        k1k2_pair = "{}-{}".format(mode_labels[K_range[i]][j],mode_labels[K_range[i+1]][jj])
                        cost = float(cost_acrossK_cons[k1k2_pair])
                        edge_w = 0.9-(cost-cost_min)/(cost_max-cost_min)*0.7
                        
                        cpatch = ConnectionPatch((.5, 0), (.5, 1), "axes fraction", "axes fraction",
                          axesA=axes[i,j], axesB=axes[i+1,jj], zorder=0, 
                          color=plt.cm.Greys(edge_w), linewidth=5, alpha=0.8)
                        axes[0,0].add_artist(cpatch)

                        coord1 = transFigure.transform(axes[i,j].transData.transform([N//2,-0.05]))
                        coord2 = transFigure.transform(axes[i+1,jj].transData.transform([N//2,1.05]))

                        text_frac = 0.2+(jj%2)*0.6 # make the text sparse when edges are dense
                        plt.text(get_center_coord(coord1[0],coord2[0],text_frac),get_center_coord(coord1[1],coord2[1],text_frac), 
                            "{:.3f}".format(cost), fontsize=18, transform=fig.transFigure, bbox=textbox_props, 
                            zorder=9, verticalalignment='center',horizontalalignment='center')
        else:
            transFigure = fig.transFigure.inverted()
            cost_max = max(cost_acrossK_cons.values())
            cost_min = min(cost_acrossK_cons.values())
            N = Q.shape[0]

            for i in range(n_row-1):
                k1k2_pair = "{}-{}".format(mode_labels[K_range[i]][0],mode_labels[K_range[i+1]][0])
                cost = float(cost_acrossK_cons[k1k2_pair])
                edge_w = 0.9-(cost-cost_min)/(cost_max-cost_min)*0.8
                
                cpatch = ConnectionPatch((.5, 0), (.5, 1), "axes fraction", "axes fraction",
                  axesA=axes[i], axesB=axes[i+1], zorder=0, 
                  color=plt.cm.Greys(edge_w), linewidth=5, alpha=0.8)
                axes[0].add_artist(cpatch)
                coord1 = transFigure.transform(axes[i].transData.transform([N//2,-0.05]))
                coord2 = transFigure.transform(axes[i+1].transData.transform([N//2,1.05]))
                text_frac = 0.5 # make the text sparse when edges are dense
                plt.text(get_center_coord(coord1[0],coord2[0],text_frac),get_center_coord(coord1[1],coord2[1],text_frac), 
                    "{:.3f}".format(cost), fontsize=18, transform=fig.transFigure, bbox=textbox_props, 
                    zorder=9, verticalalignment='center',horizontalalignment='center')


        fig.savefig(os.path.join(fig_path,"modes_aligned_multipartite_{}.png".format(cons_suffix)), bbox_inches='tight',dpi=300)
        plt.close(fig)
    

   

# get colormap for visualization
def get_random_cmap(base_cmap,K_max,seed=999):
    np.random.seed(seed)
    cmap = base_cmap(np.linspace(0, 1, K_max))
    np.random.shuffle(cmap)
    cmap = cm.colors.ListedColormap(cmap)
    return cmap

def plot_colorbar(cmap,K_max,fig_path):
    colors = cmap(np.arange(K_max))
    fig, ax = plt.subplots(figsize=(1.5*(K_max-1), 2),facecolor='white',
                           subplot_kw=dict(yticks=[])) #xticks=[], 
    plt.xticks(ticks=np.arange(0.5,K_max+0.5,1),labels=["Cls.{}".format(k) for k in np.arange(1,K_max+1)])
    ax.imshow([colors], extent=[0, K_max, 0, 1])
    # plt.title("color of each cluster")
    fig.savefig(os.path.join(fig_path,"colorbar.png"), bbox_inches='tight',dpi=300, transparent=False)
    plt.close(fig)


def plot_membership(ax,P,K_max,cmap,title):

    N = P.shape[0]
    P_aug = np.hstack((np.zeros((N,1)),P))

    K = P.shape[1]
    for k in range(K):
        ax.bar(range(N), P_aug[:,(k+1)], bottom=np.sum(P_aug[:,:(k+1)],axis=1), 
               width=1.0, edgecolor='w', linewidth=0, facecolor=cmap(k/K_max))

    ax.set_xticks([])
    ax.set_xlim([-0.5,N-0.5])
    ax.set_ylim([0,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([])
    if title:
        ax.set_ylabel("\n".join(title.split()), rotation=0, fontsize=18, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    return


def plot_membership_sorted(ax,P,K_max,cmap,title):

    N = P.shape[0]
    P_aug = np.hstack((np.zeros((N,1)),P))

    K = P.shape[1]
    for i in range(N):
      q = P[i,:]
      indices = np.argsort(-q)
      q_sorted = np.concatenate(([0],q[indices]))
      for k in range(K):
        ax.bar(i, q_sorted[(k+1)], bottom=np.sum(q_sorted[:(k+1)]), 
                  width=1.0, edgecolor='w', linewidth=0, facecolor=cmap(indices[k]/K))

    ax.set_xticks([])
    ax.set_xlim([0,N])
    ax.set_ylim([0,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([])
    # ax.set_title(title)
    if title:
        ax.set_ylabel("\n".join(title.split()), rotation=0, fontsize=18, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    return


def plot_withinK_modes(K,K_max,consQ_modes,alignment_acrossK,fig_path,cmap,fig_suffix=""):
    
    modes = range(len(consQ_modes[K]))
    fig, axes = plt.subplots(len(modes),1,figsize=(20,2*len(modes)),facecolor='white')
    if len(modes)==1:
        ax = axes
    else:
        ax = axes[0]

    m = 0
    lb1 = "K{}M{}".format(K,m+1)
    P = consQ_modes[K][m]
    plot_membership(ax,P,K_max,cmap,lb1)
    
    
    for i in range(1,len(modes)):
        m = modes[i]
        lb2 = "K{}M{}".format(K,m+1)
        # retrieve alignment pattern
        alignment = alignment_acrossK["{}-{}".format(lb1,lb2)]
        Q = consQ_modes[K][m]
        aligned_Q = alignQ_wrtP(Q,Q,alignment,merge=True)
        # plot
        ax = axes[i]
        plot_membership(ax,aligned_Q,K_max,cmap,lb2)
    
    fig.savefig(os.path.join(fig_path,"K{}_modes_{}.png".format(K,fig_suffix)), bbox_inches='tight',dpi=300)
    plt.close(fig)
    
    
def alignQ_wrtP_pattern(P,Q,idxQ2P,merge=True):
    idxQ2P = list(idxQ2P)
    if merge:        
        aligned_Q = np.zeros_like(P)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        dups = np.unique([i for i in idxQ2P if idxQ2P.count(i)>1])
        extras = list()
        dups_min = defaultdict(lambda: (float('inf'),None))

        new_pattern = [0 for _ in range(Q.shape[1])]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx not in dups:
                new_pattern[q_idx] = p_idx
                aligned_Q[:,p_idx] = Q[:,q_idx]
            else:
                diff = np.linalg.norm(Q[:,q_idx]-P[:,p_idx])
                if dups_min[p_idx][0] > diff:
                    dups_min[p_idx] = (diff,q_idx) 
        extra_cnt = P.shape[1]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx in dups:
                if q_idx==dups_min[p_idx][1]:
                    new_pattern[q_idx] = p_idx
                    aligned_Q[:,p_idx] = Q[:,q_idx]
                else:
                    new_pattern[q_idx] = extra_cnt
                    extras.append(q_idx)
                    extra_cnt += 1
        
        for ie, e in enumerate(extras):
            aligned_Q[:,P.shape[1]+ie] = Q[:,e]
            
    return aligned_Q, new_pattern


def plot_replicates(Q_list,K_range,output_path,cmap):

    K_max = max(K_range)
    num_figs = len(Q_list)
    fig_path = os.path.join(output_path,"visualization")

    fig, axes = plt.subplots(num_figs//2,2,figsize=(20,num_figs//2),facecolor='white')
    for fig_idx in range(num_figs):
        plot_membership(axes[fig_idx//2,fig_idx%2],Q_list[fig_idx],K_max,cmap,"") 

    fig.savefig(os.path.join(fig_path,"all_replicates_unaligned.png"), bbox_inches='tight',dpi=300)
    plt.close(fig)

    return


def plot_major_modes(K_range,output_path,cons_suffix,cmap):
    
    num_major_modes = len(K_range) 
    K_max = max(K_range)
    base_patterns = dict()

    fig_path = os.path.join(output_path,"visualization")
    modeQ_path = os.path.join(output_path,"modes_aligned")

    fig, axes = plt.subplots(num_major_modes,1,figsize=(15,1.5*num_major_modes),facecolor='white')
    fig_idx = 0

    for K in K_range:
        Q = np.loadtxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format("K{}M1".format(K),cons_suffix)), delimiter=' ')
        plot_membership(axes[fig_idx],Q,K_max,cmap,"K={}".format(K)) 
        fig_idx += 1

    fig.savefig(os.path.join(fig_path,"major_modes_aligned_{}.png".format(cons_suffix)), bbox_inches='tight',dpi=300)
    plt.close(fig)

    return

def plot_all_modes(K_range,mode_labels,output_path,cons_suffix,cmap):
    
    num_modes = sum([len(mode_labels[k]) for k in K_range]) 
    K_max = max(K_range)
    base_patterns = dict()

    fig_path = os.path.join(output_path,"visualization")
    modeQ_path = os.path.join(output_path,"modes_aligned")

    fig, axes = plt.subplots(num_modes,1,figsize=(15,1.5*num_modes),facecolor='white')
    fig_idx = 0

    for K in K_range:
        for i_m,m in enumerate(mode_labels[K]):
            Q = np.loadtxt(os.path.join(modeQ_path,'{}_aligned_{}.Q'.format(m,cons_suffix)), delimiter=' ')
            plot_membership(axes[fig_idx],Q,K_max,cmap,"K={}".format(K) if i_m==0 else "") 
            fig_idx += 1

    fig.savefig(os.path.join(fig_path,"all_modes_aligned_{}.png".format(cons_suffix)), bbox_inches='tight',dpi=300)
    plt.close(fig)

    return
 
    
def load_Q(Q_path,Q_files):

    Q_list = list()
    K_list = list()
       
    R = len(Q_files)
    if R==0:
        # sanity check if data directory is empty
        sys.exit("ERROR: no Q files detected. Please double check input_path and input_format.")
    for r in range(R):
        Q = np.loadtxt(os.path.join(Q_path,Q_files[r]))        
        K = Q.shape[1]
        if r==0:
            N = Q.shape[0]
        else:
            assert(Q.shape[0]==N)
            
        Q_list.append(Q)
        K_list.append(K)
        
    # reorder by consecutive K values   
    sored_idx = list(np.argsort(K_list))
    K_list = [K_list[i] for i in sored_idx]
    Q_list = [Q_list[i] for i in sored_idx]
    file_list = [Q_files[i] for i in sored_idx]
    N = Q.shape[0]

    return N, R, Q_list, K_list, file_list

def load_ind(file_path):
        
    info = pd.read_csv(os.path.join(file_path,'ind_info.txt'),delimiter=" ",index_col=0)
    # get populization sizes
    popIDs = info["Pop"].unique()
    pop_n_ind = {}
    for p in popIDs:
        pop_n_ind[p] = sum(np.array(info["Pop"]==p))
    ind2pop = info["Pop"].values
    
    return ind2pop, pop_n_ind

def load_Q_and_indinfo(input_base_path,cons_suffix,indinfo=True):
    
    input_names = [f for f in os.listdir(input_base_path) if os.path.isdir(os.path.join(input_base_path, f))]
    
    # load mode Q and ind pop info    
    N_all = list()
    R_all = list()
    Q_all = list()
    K_all = list()
    if indinfo:
        ind2pop_all = list()
        popNind_all = list()

    for input_name in input_names:
        Q_path = os.path.join(input_base_path, input_name,"modes_aligned")
        Q_files = [i for i in os.listdir(Q_path) if i.endswith('{}.Q'.format(cons_suffix))]
        N, R, Q_list, K_list, file_list = load_Q(Q_path,Q_files)
        
        N_all.append(N)
        R_all.append(R)
        Q_all.append(Q_list)
        K_all.append(K_list)

        if indinfo:
            ind2pop, pop_n_ind = load_ind(os.path.join(input_base_path,input_name,"input"))
            ind2pop_all.append(ind2pop)
            popNind_all.append(pop_n_ind)
    if indinfo:    
        return input_names,N_all, R_all, Q_all, K_all, ind2pop_all, popNind_all
    else:
        return input_names,N_all, R_all, Q_all, K_all


def plot_membership_with_pop(ax,P,ind2pop,K_max,cmap,title):

    N = P.shape[0]
    P_aug = np.hstack((np.zeros((N,1)),P))

    K = P.shape[1]
    for k in range(K):
        ax.bar(range(N), P_aug[:,(k+1)], bottom=np.sum(P_aug[:,0:(k+1)],axis=1), 
               width=1.0, edgecolor='w', linewidth=0, facecolor=cmap(k/K_max))
    # pop_split = np.where(np.diff(ind2pop)==1)[0]
    # ax.vlines(pop_split+0.5, 0, 1, colors='gray', linewidth=0.6, alpha=0.3)
    
    ax.set_xticks([])
    ax.set_xlim([0,N])
    ax.set_ylim([0,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([])
    # ax.set_title(title)
    if title:
        ax.set_ylabel("\n".join(title.split()), rotation=0, fontsize=18, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    return

def align_popwise_membership(input_names, Q_all, K_all, ind2pop_all, output_path, cmap):

    K_range = np.sort(list(set().union(*[set(K) for K in K_all])))
    K_max = np.max(K_range)
    
    num_input = len(input_names)
    common_pops = list(set.intersection(*[set(inds) for inds in ind2pop_all]))
    n_pop = len(common_pops)
    
    # align modes using population-wise memberships
    cost_res = dict()
    align_res = dict()
    
    with open(os.path.join(output_path,"popwise_inputs.txt"),"w") as f:
        for i,input in enumerate(input_names):
            f.write("{},{}\n".format(i+1,input))
    
    with open(os.path.join(output_path,"popwise_alignment.txt"),"w") as f:
        f.write("K,input1_index,input2_index,optimal_obj,alignment\n")
        
        for K in K_range[::-1]:
            # align each K (to input 0)
            Qpop_list = list()
            pop_size_list = list()
            for q in range(num_input):
                Q = Q_all[q][np.where(np.array(K_all[q])==K)[0][0]]
                ind2pop = ind2pop_all[q]
                # summarize Q by pop 
                Qpop = np.zeros((n_pop,K))
                pop_size = np.zeros((n_pop,))
                for i_pop in range(n_pop):
                    Qpop[i_pop,:] = np.mean(Q[ind2pop==common_pops[i_pop],:],axis=0)
                    pop_size[i_pop] = np.sum(ind2pop==common_pops[i_pop])
                Qpop_list.append(Qpop)
                pop_size_list.append(pop_size)
            
            cost_res[K] = [[np.nan for _ in input_names] for _ in input_names]  
            align_res[K] = [[None for _ in input_names] for _ in input_names]    
            
            for i_p in range(num_input):
                P = Qpop_list[i_p]
                popsize_p = pop_size_list[i_p]
                for i_q in range(num_input):
                    Q = Qpop_list[i_q]
                    popsize_q = pop_size_list[i_q]
                    weight = (popsize_p+popsize_q)/np.sum(popsize_p+popsize_q)
                                
                    opt_obj, idxQ2P = align_ILP_weighted(P, Q, weight)
                    
                    f.write("{},{},{},{},{}\n".format(K,i_p+1,i_q+1,opt_obj," ".join([str(id+1) for id in idxQ2P])))
                    
                    cost_res[K][i_p][i_q] = opt_obj
                    align_res[K][i_p][i_q] = idxQ2P
    
    # plot alignment
    
    K_maxcnt = defaultdict(lambda: float('-inf'))
    for d in [Counter(K) for K in K_all]:
        for k, v in d.items():
            K_maxcnt[k] = max(K_maxcnt[k], v)
    K_cnt_list = [K_maxcnt[k] for k in K_range]
    tot_K_cnt = np.sum(K_cnt_list)
    sep_K_cnt = np.cumsum(K_cnt_list)
    sep_K_cnt = np.insert(sep_K_cnt, 0, 0)
    
    fig, axes = plt.subplots(tot_K_cnt,num_input,figsize=(30,2*(tot_K_cnt)),facecolor='white')
    
    for i_K,K in enumerate(K_range):
        i_p = 0
        k_idx = 0
        for idx in np.where(np.array(K_all[i_p])==K)[0]:
            P = Q_all[i_p][idx]
            ax = axes[k_idx+sep_K_cnt[i_K],i_p]
            plot_membership_with_pop(ax,P,ind2pop_all[i_p],K_max,cmap,"K={}".format(K) if k_idx==0 else "")
            k_idx += 1
        while k_idx<K_maxcnt[K]:
            ax = axes[k_idx+sep_K_cnt[i_K],i_p]
            ax.axis('off')
            k_idx += 1
        
        for i_q in range(1,num_input):
            k_idx = 0
            for idx in np.where(np.array(K_all[i_q])==K)[0]:
                Q = Q_all[i_q][idx]
                idxQ2P = align_res[K][i_p][i_q]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
                ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                plot_membership_with_pop(ax,aligned_Q,ind2pop_all[i_q],K_max,cmap,"")
                k_idx += 1
            while k_idx<K_maxcnt[K]:
                ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                ax.axis('off')
                k_idx += 1
        
    for i_p in range(num_input):
        ax = axes[0,i_p]
        ax.set_title(input_names[i_p], fontsize=20)
        
    fig.savefig(os.path.join(output_path,"aligned_all.png"), bbox_inches='tight',dpi=300)
    plt.close(fig)


def align_multiple_model(input_names, Q_all, K_all, output_path, cmap):

    K_range = np.sort(list(set().union(*[set(K) for K in K_all])))
    K_max = np.max(K_range)
    
    num_input = len(input_names)
    
    cost_res = dict()
    align_res = dict()
    
    with open(os.path.join(output_path,"multiple_inputs.txt"),"w") as f:
        for i,input in enumerate(input_names):
            f.write("{},{}\n".format(i+1,input))
    
    with open(os.path.join(output_path,"multiple_alignment.txt"),"w") as f:
        f.write("K,input1_index,input2_index,optimal_obj,alignment\n")
        
        for K in K_range[::-1]:
            Qrep_list = list()
            # align each K (to input 0)
            for q in range(num_input):
                Q = Q_all[q][np.where(np.array(K_all[q])==K)[0][0]]
                Qrep_list.append(Q)
            
            cost_res[K] = [[np.nan for _ in input_names] for _ in input_names]  
            align_res[K] = [[None for _ in input_names] for _ in input_names]    
            
            for i_p in range(num_input):
                P = Qrep_list[i_p]

                for i_q in range(num_input):
                    Q = Qrep_list[i_q]
                                
                    opt_obj, idxQ2P = align_ILP(P, Q)
                    
                    f.write("{},{},{},{},{}\n".format(K,i_p+1,i_q+1,opt_obj," ".join([str(id+1) for id in idxQ2P])))
                    
                    cost_res[K][i_p][i_q] = opt_obj
                    align_res[K][i_p][i_q] = idxQ2P
                
    
    # plot alignment
    K_maxcnt = defaultdict(lambda: float('-inf'))
    for d in [Counter(K) for K in K_all]:
        for k, v in d.items():
            K_maxcnt[k] = max(K_maxcnt[k], v)
    K_cnt_list = [K_maxcnt[k] for k in K_range]
    tot_K_cnt = np.sum(K_cnt_list)
    sep_K_cnt = np.cumsum(K_cnt_list)
    sep_K_cnt = np.insert(sep_K_cnt, 0, 0)
    
    fig, axes = plt.subplots(tot_K_cnt,num_input,figsize=(30,2*(tot_K_cnt)),facecolor='white')
    
    for i_K,K in enumerate(K_range):
        i_p = 0
        k_idx = 0
        for idx in np.where(np.array(K_all[i_p])==K)[0]:
            P = Q_all[i_p][idx]
            ax = axes[k_idx+sep_K_cnt[i_K],i_p]
            plot_membership(ax,P,K_max,cmap,"K={}".format(K) if k_idx==0 else "")
            k_idx += 1
        while k_idx<K_maxcnt[K]:
            ax = axes[k_idx+sep_K_cnt[i_K],i_p]
            ax.axis('off')
            k_idx += 1
        
        for i_q in range(1,num_input):
            k_idx = 0
            for idx in np.where(np.array(K_all[i_q])==K)[0]:
                Q = Q_all[i_q][idx]
                idxQ2P = align_res[K][i_p][i_q]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
                ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                plot_membership(ax,aligned_Q,K_max,cmap,"")
                k_idx += 1
            while k_idx<K_maxcnt[K]:
                ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                ax.axis('off')
                k_idx += 1
        
    for i_p in range(num_input):
        ax = axes[0,i_p]
        ax.set_title(input_names[i_p], fontsize=20)
        
    fig.savefig(os.path.join(output_path,"aligned_all.png"), bbox_inches='tight',dpi=300)
    plt.close(fig)