import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import time
from itertools import product,combinations_with_replacement, combinations
from collections import defaultdict
from scipy import special

import networkx as nx
from TracyWidom import TracyWidom
import community.community_louvain as community_louvain
from scipy.spatial.distance import cdist
import cvxpy as cp

from clumppling.func_utils import *

#%% Main Functions

def align_ILP(P,Q):
    
    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    N_ind = Q.shape[0] #individuals
    assert(K_Q>=K_P)
    
    # compute distance
    D = cdist(Q.T, P.T, metric='sqeuclidean')/N_ind # seuclidean, euclidean : 
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
            # for i,idx in enumerate(old_idx):
            idxQ2P[old_idx] = idxQ2P_comb[np.arange(len(old_idx))]

            idxQ2P[pair[0]] = idxQ2P_comb[new_idx]
            idxQ2P[pair[1]] = idxQ2P_comb[new_idx]
    idxQ2P = idxQ2P.astype(int)


    return opt_obj, idxQ2P

# def align_ILP_weighted(P,Q,weight):
    
#     K_Q = Q.shape[1] # #rows
#     K_P = P.shape[1] # #columns
#     N_ind = Q.shape[0] #individuals
#     assert(K_Q>=K_P)

#     # compute distance
#     dist_elem = (np.expand_dims(Q, axis=2)-np.expand_dims(P, axis=1))**2
    
#     w_dist = dist_elem*np.expand_dims(weight,(1,2))
#     D = np.sum(w_dist,axis=0)
#     d = D.flatten()
    
#     # constraints
#     # sum of row = 1
#     A = np.zeros((K_Q,d.shape[0]))
#     for i in range(K_Q):
#         A[i,i*K_P:(i+1)*K_P] = np.ones((1,K_P))
#     A_ones = np.ones((K_Q))
    
#     # sum of column >= 1
#     B = np.zeros((K_P,d.shape[0]))
#     for i in range(K_P):
#         B[i,np.arange(i,d.shape[0],K_P)] = np.ones((1,K_Q))
#     B_ones = np.ones((K_P))
 
#     # ILP
#     w = cp.Variable(d.shape[0], boolean = True)
#     constraints = [A @ w == A_ones, B @ w >= B_ones]
#     obj = d.T @ w
#     problem = cp.Problem(cp.Minimize(obj), constraints)
#     problem.solve(solver=cp.GLPK_MI) #
    
#     # optimal solution
#     opt_obj = problem.value
#     opt_w = problem.variables()[0].value.astype(int)
#     opt_W = np.reshape(opt_w,D.shape)
#     # print("Optimal value: {}".format(opt_obj))
#     # print("Solution:\n {}".format(opt_W))
    
#     # matching index Q_idx to P_idx
#     idxQ2P = np.where(opt_W==1)[1]
    
    
#     return opt_obj, idxQ2P#, idxP2Q


def alignQ_wrtP(P,Q,idxQ2P,merge=True):
    # K1, K2 = P.shape[1], Q.shape[1]
    
    if merge:        
        aligned_Q = np.zeros_like(P)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        dups = np.unique([i for i in idxQ2P if idxQ2P.count(i)>1])

        # extras = list()
        extra_cnt = 0
        dups_min = defaultdict(lambda: (float('inf'),None))

        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx not in dups:
                aligned_Q[:,p_idx] = Q[:,q_idx]
            else:
                diff = np.linalg.norm(Q[:,q_idx]-P[:,p_idx])
                if dups_min[p_idx][0] > diff:
                    dups_min[p_idx] = (diff,q_idx) 
        
        P_dim2 = P.shape[1]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx in dups:
                if q_idx==dups_min[p_idx][1]:
                    aligned_Q[:,p_idx] = Q[:,q_idx]
                else:
                    # extras.append(q_idx)
                    aligned_Q[:,P_dim2+extra_cnt] = Q[:,q_idx]
                    extra_cnt += 1
            
    return aligned_Q


def align_ILP_withinK(ILP_withinK_filename,output_path,Q_list,K_range,k2ids):

    f = open(os.path.join(output_path,ILP_withinK_filename),"w")

    for K in K_range:
        ids = k2ids[K]
        n_ids = len(ids)
             
        for i in range(n_ids-1):
            P = Q_list[ids[i]]
            
            for j in range(i+1,n_ids):
    
                Q = Q_list[ids[j]]
                opt_obj, idxQ2P = align_ILP(P, Q)
                f.write("{} {} {} {}\n".format(K,ids[i],ids[j],opt_obj))
                f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))  

    f.close()
    

def load_ILP_withinK(Q_list,ILP_withinK_filename,output_path,K_range,k2ids,idx2idxinK,get_cost=True):

    f = open(os.path.join(output_path,ILP_withinK_filename),"r")
    alignments = f.readlines()
    f.close()
    
    # load alignments
    cost_ILP_res = dict()
    align_ILP_res = dict()
    
    for K in K_range:
        ids = k2ids[K]
        n_ids = len(ids)
         
        # align_ILP_res[K] = [[None for _ in ids] for _ in ids]
        # cost_ILP_res[K] = [[np.nan for _ in ids] for _ in ids]
        align_ILP_res[K] = np.zeros((n_ids,n_ids,K)).astype(int)
        cost_ILP_res[K] = np.zeros((n_ids,n_ids))
        
    for l in range(len(alignments)//2):
        l1 = alignments[2*l].split()
        l2 = alignments[2*l+1].split()
        
        # cost
        K = int(l1[0])
        ids = k2ids[K]
        i,j = idx2idxinK[int(l1[1])], idx2idxinK[int(l1[2])]
        
        # alignment
        # l2 = [int(idx) for idx in l2]
        # align_ILP_res[K][i][j] = [l2.index(idx) for idx in range(K)] # P2Q
        # align_ILP_res[K][j][i] = l2 # Q2P

        l2 = [int(idx) for idx in l2]
        align_ILP_res[K][i][j] = np.array([l2.index(idx) for idx in range(K)]) # P2Q
        align_ILP_res[K][j][i] = np.array(l2) # Q2P
        
        if get_cost:  
            # alignment cost
            P = Q_list[int(l1[1])]
            Q = Q_list[int(l1[2])]
            aligned_idxQ2P = l2
            aligned_Q = np.zeros_like(P)
            for q_idx in range(Q.shape[1]):
                aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
            cost_actual = cost_membership(P,aligned_Q,P.shape[0])
            cost_ILP_res[K][i][j] = cost_actual 
            cost_ILP_res[K][j][i] = cost_actual 
        
    if get_cost:
        return align_ILP_res, cost_ILP_res
    else:
        return align_ILP_res

def cd_default(G):
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

    resolution = 1
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
    """

    # Please comment out the following line and customize your community detection method here.
    partition_map = {i:0 for i in range(G.number_of_nodes())} 
    
    return partition_map

def detect_modes(cost_ILP_res,K_range,default_cd,cd_mod_thre, save_path=None):

    modes_allK_list = dict()
    msg = ""
            
    for K in K_range:
    
        # normalize cost matrix and obtain adjacency matrix
        # cost_mat = np.array(cost_ILP_res[K])
        cost_mat = cost_ILP_res[K]
        n_nodes = cost_mat.shape[0]
        if np.nanmean(cost_mat)<1e-6:
            # adj_mat = np.zeros(cost_mat.shape)
            partition_map = {i:0 for i in range(n_nodes)}         
        else:
            adj_mat = get_adj_mat(cost_mat)

            # test for community structure
            has_comm_struc = test_comm_struc(adj_mat, alpha = 0.01)
            if not has_comm_struc:
                partition_map = {i:0 for i in range(n_nodes)} 
                msg = "K={}: no significant community structure -> set to single mode".format(K)

            else: # having community structure

                # create graph
                G = nx.from_numpy_matrix(adj_mat)   
                G.remove_edges_from(nx.selfloop_edges(G))
                pos = nx.spring_layout(G)
                    
                # adj_mat_allK_list[K] = adj_mat
                     
                ########################################
                # community detection to find modes
                ########################################
                if np.nanmean(cost_mat)<1e-6:
                    partition_map = {i:0 for i in range(n_nodes)} 
                else:
                    if default_cd:
                        partition_map = cd_default(G)
                        # print("Running default Louvain for mode detection.")
                    else:
                        partition_map = cd_custom(G)
                        # print("Running customized community detection.")
                    if cd_mod_thre!=-1:
                        mod_res = community_louvain.modularity(partition_map, G)
                        if mod_res<cd_mod_thre: # quality of community detection is low --> no community structure
                            msg = "K={}: low community detection quality (below modularity threshold) -> set to single mode".format(K)
                            partition_map = {i:0 for i in range(G.number_of_nodes())} 
                
        ########################################
        # separate modes
        ########################################
        
        modes = defaultdict(list)
        for r in partition_map.keys():
            modes[partition_map[r]].append(r)
        
        modes_allK_list[K] = modes
        
    return modes_allK_list, msg
    
        
        # elif method_name=="mcl":
            
        #     inflation = method_para if method_para else 2
        #     #     coms = algorithms.markov_clustering(G,inflation=inflation,expansion=2)
        #     #     partition = coms.communities
        #     result = mc.run_mcl(nx.to_numpy_matrix(G),inflation=inflation)     # run MCL with default parameters
        #     partition = mc.get_clusters(result)
        #     partition = [set(a) for a in partition]
        #     partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}        
        
        # elif method_name=="agdl":
        #     kc = method_para if method_para else 3
        #     # choose number of cluster that yields largest modularity
        #     max_mod = -np.inf
        #     partition = None
        #     for ncl in range(1,G.number_of_nodes()//3+1):
        #         coms = algorithms.agdl(G,number_communities=ncl,kc=kc)
        #         mod = evaluation.newman_girvan_modularity(G,coms).score
        #         # print(ncl,mod)
        #         if mod>max_mod:
        #             max_mod = mod
        #             partition = coms.communities
        #     partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
        
        # elif method_name=="mcode":
        #     coms = algorithms.mcode(G,"weight")
        #     partition = coms.communities
        #     partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
        
        # else:
        #     func = getattr(algorithms, 'belief')
        #     coms = func(G)
        #     partition = coms.communities
        #     partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
            
        # else: # no method
        #     print("K={}: no algorithm provided, set to single mode".format(K))
        #     partition_map = {i:0 for i in range(G.number_of_nodes())}       
    
    

def extract_modes(Q_list,modes_allK_list,align_ILP_res,cost_ILP_res,K_range,N_ind,k2ids,save_modes=False, save_path=None,ILP_modes_filename=None):
 
    rep_modes = defaultdict(list)
    repQ_modes = dict()
    meanQ_modes = dict()
    average_stats = dict()
    
    ########################################
    # get results of alignment for each K
    ########################################
    
    for K in K_range:
    
        modes = modes_allK_list[K]
        
        # cost_mat = np.array(cost_ILP_res[K])
        cost_mat = cost_ILP_res[K]
        adj_mat = get_adj_mat(cost_mat)
    
        mds = range(len(modes.keys()))
        mds_indices = list()
        
        meanQ_modes[K] = dict()
        repQ_modes[K] = dict()
        average_stats[K] = dict()

        for mode_idx in mds:
            # print("mode {}".format(mode_idx))
    
            ########################################
            # find the concensus membership in mode
            ########################################
            all_indices = modes[mode_idx]
            if len(all_indices)==1:
                idx = all_indices[0]
                rep_modes[K].append(idx)
                average_stats[K][mode_idx] = {"cost":0,"Hprime":1} 
                minC_rep_idx = idx+k2ids[K][0]
                P = Q_list[minC_rep_idx]
                repQ_modes[K][mode_idx] = P
                meanQ_modes[K][mode_idx] = P
                
            else:
            
                sub_adj = adj_mat[np.array(all_indices),:][:,np.array(all_indices)]
                
                # find the one with max row sum (min cost)
                minC_idx = np.argmax(sub_adj.sum(axis=0))
                minC_self_idx = all_indices[minC_idx]
                minC_rep_idx = all_indices[minC_idx]+k2ids[K][0]
                rep_modes[K].append(minC_rep_idx)
                mds_indices.append(all_indices[minC_idx])
                # other non-representative ones
                other_self_indices = [m for im,m in enumerate(all_indices) if im!=minC_idx]
                other_rep_indices = [m+k2ids[K][0] for im,m in enumerate(all_indices) if im!=minC_idx]
                
                ########################################
                # compute the average cost/similarity
                ########################################
                
                comm_cost_mat = cost_mat[np.array(all_indices),:][:,np.array(all_indices)]
                Gprime_mat = C2Gprime(comm_cost_mat)
                avg_cost = avg_tril(comm_cost_mat)
                h_prime = avg_tril(Gprime_mat)
                average_stats[K][mode_idx] = {"cost":avg_cost,"Hprime":h_prime}     
            
                ########################################
                # get mean Q over all in the mode
                ########################################
                P = Q_list[minC_rep_idx]
                repQ_modes[K][mode_idx] = P
                
                Q_sum = np.zeros_like(P)
                Q_sum += Q_list[minC_rep_idx]
                # Q_aligned_list = list()
                for i,r in enumerate(other_rep_indices):
                    Q = Q_list[r]
                    aligned_idxQ2P = align_ILP_res[K][other_self_indices[i]][minC_self_idx]
                    aligned_Q = np.zeros_like(P)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
                    Q_sum += aligned_Q

                meanQ_modes[K][mode_idx] = Q_sum/len(all_indices)

    if save_modes:
        # write the by-mode alignment into file
        file_name = os.path.join(save_path,ILP_modes_filename)
        write_modes_to_file(file_name,K_range,N_ind,modes_allK_list,meanQ_modes)

    return rep_modes,repQ_modes,meanQ_modes,average_stats   


def align_leader_clustering(cost_thre,Q_list,K_range,N_ind,k2ids,save_modes=False, save_path=None,ILP_modes_filename=None):

    modes_allK_list = dict()
    rep_modes = defaultdict(list)
    repQ_modes = dict()
    meanQ_modes = dict()
    average_stats = dict()
    align_ILP_res = dict()


    for K in K_range:
        
        repQ_modes[K] = dict()
        meanQ_modes[K] = dict()
        average_stats[K] = dict()
        align_in_K = defaultdict(dict)
        
        ids = k2ids[K]
        n_ids = len(ids)

        leaders = list()
        leader_followers = dict()
        lc_cost = dict()
        lc_aligned = dict()
        
        l = 0
        leaders.append(l)
        leader_followers[l] = list()
        lc_aligned[l] = list()
        lc_cost[l] = list()
        
        
        for r in np.arange(1,n_ids):
            Q = Q_list[ids[r]]
            
            perm_to_leaders = list()
            cost_to_leaders = list()
            
            for l in leaders:
                P = Q_list[ids[l]]
                opt_obj, idxQ2P = align_ILP(P, Q)
                cost_to_leaders.append(opt_obj)
                perm_to_leaders.append(idxQ2P)
            
            min_cost_idx = np.argmin(cost_to_leaders)
            min_cost = cost_to_leaders[min_cost_idx]

            if min_cost>(cost_thre): #*1.2**(K-K_range[0])
                leaders.append(r)
                leader_followers[r] = list()
                lc_aligned[r] = list()
                lc_cost[r] = list()
            else:
                l = leaders[min_cost_idx]
                leader_followers[l].append(r)
                lc_aligned[l].append(perm_to_leaders[min_cost_idx])
                lc_cost[l].append(min_cost)

        # add to modes
        modes = defaultdict(list)
        for i_l,l in enumerate(leaders):
            modes[i_l].append(l)
            rep_modes[K].append(ids[l])
            P = Q_list[ids[l]]
            repQ_modes[K][i_l] = P
            
            
            Q_sum = np.zeros_like(P)
            Q_sum += P
            
            ########################################
            # compute the average cost/similarity
            ########################################
            if len(leader_followers[l])==0:
                average_stats[K][i_l] = {"cost":0,"Hprime":1} 
                
            else:
                avg_cost = np.mean(lc_cost[l])
                h_prime = np.mean(C2Gprime(lc_cost[l]))
                average_stats[K][i_l] = {"cost":avg_cost,"Hprime":h_prime}     
            
                ############################################
                # get modes and mean Q over all in the mode
                ############################################
        
                for i_fol, follower in enumerate(leader_followers[l]):
                    modes[i_l].append(follower)
                    
                    Q = Q_list[ids[follower]]
                    aligned_idxQ2P = lc_aligned[l][i_fol]
                    align_in_K[follower][l] = aligned_idxQ2P
                    
                    aligned_Q = np.zeros_like(P)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
                    Q_sum += aligned_Q
                
            meanQ_modes[K][i_l] = Q_sum/(len(leader_followers[l])+1)
        
        align_ILP_res[K] = align_in_K
        modes_allK_list[K] = modes  
        
    if save_modes:
        # write the by-mode alignment into file
        file_name = os.path.join(save_path,ILP_modes_filename)
        write_modes_to_file(file_name,K_range,N_ind,modes_allK_list,meanQ_modes)
        
    
    
    return modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats 

def align_leader_clustering_adaptive(cost_thre_base,Q_list,K_range,N_ind,k2ids,ind2pop, pop_n_ind,save_modes=False, save_path=None,ILP_modes_filename=None):

    modes_allK_list = dict()
    rep_modes = defaultdict(list)
    repQ_modes = dict()
    meanQ_modes = dict()
    average_stats = dict()
    align_ILP_res = dict()
    
    theo_cost_list = list()
    for Q in Q_list:
        theo_cost_list.append(get_theoretical_cost(Q, pop_n_ind, ind2pop))
    theo_cost_list = np.array(theo_cost_list)

    for K in K_range:
        
        repQ_modes[K] = dict()
        meanQ_modes[K] = dict()
        average_stats[K] = dict()
        align_in_K = defaultdict(dict)
        
        ids = k2ids[K]
        n_ids = len(ids)
        
        ct = np.mean(theo_cost_list[ids])

        leaders = list()
        leader_followers = dict()
        lc_cost = dict()
        lc_aligned = dict()
        # lc_theo_cost = dict()
        
        l = 0
        leaders.append(l)
        leader_followers[l] = list()
        lc_aligned[l] = list()
        lc_cost[l] = list()
        # lc_theo_cost[l] = [theo_cost_list[0]]
        
        for r in np.arange(1,n_ids):
            Q = Q_list[ids[r]]
            
            perm_to_leaders = list()
            cost_to_leaders = list()
            
            for l in leaders:
                P = Q_list[ids[l]]
                opt_obj, idxQ2P = align_ILP(P, Q)
                cost_to_leaders.append(opt_obj)
                perm_to_leaders.append(idxQ2P)
            
            min_cost_idx = np.argmin(cost_to_leaders)
            min_cost = cost_to_leaders[min_cost_idx]

            cost_thre = 5*ct #cost_thre_base+10*np.mean([theo_cost_list[ids[r]],theo_cost_list[ids[leaders[min_cost_idx]]]]) #lc_theo_cost[leaders[min_cost_idx]]
            # print(cost_thre_base,cost_thre,min_cost)
            
            if min_cost>(cost_thre): #*1.2**(K-K_range[0])
                leaders.append(r)
                leader_followers[r] = list()
                lc_aligned[r] = list()
                lc_cost[r] = list()
                # lc_theo_cost[r] = [theo_cost_list[ids[r]]]
                
            else:
                l = leaders[min_cost_idx]
                leader_followers[l].append(r)
                lc_aligned[l].append(perm_to_leaders[min_cost_idx])
                lc_cost[l].append(min_cost)
                # lc_theo_cost[l].append(theo_cost_list[ids[r]])
                
        # add to modes
        modes = defaultdict(list)
        for i_l,l in enumerate(leaders):
            modes[i_l].append(l)
            rep_modes[K].append(ids[l])
            P = Q_list[ids[l]]
            repQ_modes[K][i_l] = P
            
            
            Q_sum = np.zeros_like(P)
            Q_sum += P
            
            ########################################
            # compute the average cost/similarity
            ########################################
            if len(leader_followers[l])==0:
                average_stats[K][i_l] = {"cost":0,"Hprime":1} 
                
            else:
                avg_cost = np.mean(lc_cost[l])
                h_prime = np.mean(C2Gprime(lc_cost[l]))
                average_stats[K][i_l] = {"cost":avg_cost,"Hprime":h_prime}     
            
                ############################################
                # get modes and mean Q over all in the mode
                ############################################
        
                for i_fol, follower in enumerate(leader_followers[l]):
                    modes[i_l].append(follower)
                    
                    Q = Q_list[ids[follower]]
                    aligned_idxQ2P = lc_aligned[l][i_fol]
                    align_in_K[follower][l] = aligned_idxQ2P
                    
                    aligned_Q = np.zeros_like(P)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
                    Q_sum += aligned_Q
                
            meanQ_modes[K][i_l] = Q_sum/(len(leader_followers[l])+1)
        
        align_ILP_res[K] = align_in_K
        modes_allK_list[K] = modes  
        
    if save_modes:
        # write the by-mode alignment into file
        file_name = os.path.join(save_path,ILP_modes_filename)
        write_modes_to_file(file_name,K_range,N_ind,modes_allK_list,meanQ_modes)
            
    return modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats 


def align_ILP_modes_acrossK(consensusQ_modes,K_range,N,save_path,ILP_acrossK_filename,enum_comb_k=True,ind2pop=None):
    
    K_range_sorted = sorted(K_range,reverse=True)
    K_comb = list([(K_range_sorted[i],K_range_sorted[i+1]) for i in range(len(K_range_sorted)-1)])
    K_comb.extend([(K_range_sorted[i],K_range_sorted[i]) for i in range(len(K_range_sorted))])
    
    perm_ILP_acrossK_Q2P = dict()
    # perm_ILP_acrossK_P2Q = dict()
    cost_acrossK_ILP = dict()
    
    f = open(os.path.join(save_path,ILP_acrossK_filename),"w")
    
    best_alignments = dict()
    
    for K1,K2 in K_comb:
        
        best_alignment_idx = None
        best_alignment_obj = np.inf

        # labels = list(product(["K{}m{}".format(K1,im) for im,m in enumerate(consensusQ_modes[K1])],\
        #                       ["K{}m{}".format(K2,im) for im,m in enumerate(consensusQ_modes[K2])]))
        rijs = list(product(range(len(consensusQ_modes[K1])),range(len(consensusQ_modes[K2]))))
        
        for i,(ri,rj) in enumerate(rijs):
            
            Q = consensusQ_modes[K1][ri] # Q is the one with more clusters
            P = consensusQ_modes[K2][rj]

            if (K1==K2 and rj==ri):
                # identical alignment (dummy)
                opt_obj = 0
                idxQ2P = np.arange(K1)
            elif enum_comb_k and (K1-K2)==1:
                opt_obj, idxQ2P = align_ILP_diff1(P, Q)
            else:
                opt_obj, idxQ2P = align_ILP(P, Q)

            if best_alignment_obj>opt_obj:
                best_alignment_obj = opt_obj
                best_alignment_idx = i
            f.write("{}#{}-{}#{}:{}\n".format(K2,rj,K1,ri,opt_obj))
            f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
            # f.write("{}\n".format(" ".join([str(id) for id in idxP2Q])))
    
            cost_acrossK_ILP["{}#{}-{}#{}".format(K2,rj,K1,ri)] = opt_obj
            perm_ILP_acrossK_Q2P["{}#{}-{}#{}".format(K2,rj,K1,ri)] = idxQ2P
            # perm_ILP_acrossK_P2Q[labels[i]] = idxP2Q
        
        best_alignments[(K1,K2)] = rijs[best_alignment_idx]
    f.close()
    
    # best alignments
    best_ILP_acrossK = list()
    f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best."+ILP_acrossK_filename.split(".")[1]),"w")
    stats = ['cost', 'Hprime']
    f.write('best_pair\t{}\n'.format("\t".join(stats)))
    for i_K1 in range(len(K_range_sorted)-1):
        K1 = K_range_sorted[i_K1]
        K2 = K_range_sorted[i_K1+1]
        ri,rj = best_alignments[(K1,K2)]
        bali = "{}#{}-{}#{}".format(K2,rj,K1,ri)
        Q = consensusQ_modes[K1][ri]
        P = consensusQ_modes[K2][rj]
        aligned_Q = alignQ_wrtP(P,Q,perm_ILP_acrossK_Q2P["{}#{}-{}#{}".format(K2,rj,K1,ri)],merge=True)
        cost_bali = cost_membership(P,aligned_Q,P.shape[0])
        best_ILP_acrossK.append(bali)
        f.write("{}\t{}\t{}\n".format(bali,cost_bali,C2Gprime(cost_bali)))
    f.close()
    
    return perm_ILP_acrossK_Q2P, cost_acrossK_ILP, best_ILP_acrossK


# def load_ILP_acrossK(save_path,ILP_acrossK_filename):

#     # load alignments
#     f = open(os.path.join(save_path,ILP_acrossK_filename),"r")
#     alignments_acrossK = f.readlines()
#     f.close()
    
#     cost_acrossK_ILP_res = dict()
#     # align_acrossK_ILP_res_P2Q = dict()
#     align_acrossK_ILP_res_Q2P = dict() # Q has more clusters
    
#     n_pairs = (len(alignments_acrossK)-1)//2
    
#     for l in range(n_pairs):
#         l1 = alignments_acrossK[2*l].split(":")
#         l2 = alignments_acrossK[2*l+1].split()
#         # l3 = alignments_acrossK[3*l+2].split()
        
#         # alignment
#         l2 = [int(idx) for idx in l2]
#         # l3 = [int(idx) for idx in l3]
    
#         cost_acrossK_ILP_res[l1[0]] = l1[1]
#         # align_acrossK_ILP_res_P2Q[l1[0]] = l3 # P2Q
#         align_acrossK_ILP_res_Q2P[l1[0]] = l2 # Q2P
    
#     # best alignments
#     f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best."+ILP_acrossK_filename.split(".")[1]),"r")
#     best_ILP_acrossK = f.readlines()
#     f.close()
        
#     return align_acrossK_ILP_res_Q2P, cost_acrossK_ILP_res, best_ILP_acrossK


def report_stats(modes_allK_list,average_stats,K_range,k2ids,save_path):

    stats = ['cost', 'Hprime']
    weighted_stats = list()
    for stat in stats:
        weighted_stats.append(defaultdict(float))
    
    f = open(os.path.join(save_path,"stats.txt"),"w")
    f.write('K\tmode\tsize\t{}\n'.format("\t".join(stats)))
    for K in K_range:
        for mk in average_stats[K]:
            s = len(modes_allK_list[K][mk])
            for i_stat, stat in enumerate(stats):
                weighted_stats[i_stat][K] += s/len(k2ids[K])*average_stats[K][mk][stat]
            f.write('{}\t{}\t{}\t{}\n'.format(K,mk,s, "\t".join([str(average_stats[K][mk][stat]) for stat in stats])))
    for K in K_range:
        f.write('{}\t{}\t{}\t{}\n'.format(K,'weighted_avg',len(k2ids[K]), "\t".join([str(ws[K]) for ws in weighted_stats])))
    f.close()

    return weighted_stats


def test_comm_struc(W, alpha = 0.05):
    
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