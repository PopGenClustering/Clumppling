import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import time
from itertools import product,combinations_with_replacement
from collections import defaultdict
from scipy import special

# from cdlib import algorithms, evaluation
import community as community_louvain
import networkx as nx
# import markov_clustering as mc # MCL for mode detection

# import networkx.algorithms.community as nx_comm
# import markov_clustering as mc # MCL for mode detection

# import community.community_louvain

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
    D = cdist(Q.T, P.T, metric='euclidean')**2/N_ind # seuclidean, euclidean : 
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
    problem.solve() #solver=cp.GLPK_MI
    
    # optimal solution
    opt_obj = problem.value
    opt_w = problem.variables()[0].value.astype(int)
    opt_W = np.reshape(opt_w,D.shape)
    # print("Optimal value: {}".format(opt_obj))
    # print("Solution:\n {}".format(opt_W))
    
    # matching index Q_idx to P_idx
    idxQ2P = np.where(opt_W==1)[1]
    # idxP2Q = np.where(opt_W.T==1)[1]
    return opt_obj, idxQ2P#, idxP2Q


# def align_IQP(P,Q):
    
#     K_Q = Q.shape[1] # #rows
#     K_P = P.shape[1] # #columns
#     N_ind = Q.shape[0] #individuals
#     assert(K_Q>=K_P)
    
#     W = cp.Variable((K_Q,K_P), boolean = True)
    
#     obj = cp.sum_squares(Q @ W - P)/N_ind
    
#     constraints = [W @ np.ones((K_P,1)) == np.ones((K_Q,1)), np.ones((1,K_Q)) @ W >= np.ones((1,K_P))]
    
#     problem = cp.Problem(cp.Minimize(obj), constraints)
#     problem.solve() #solver='SCIP'
    
#     # optimal solution
#     opt_obj = problem.value
#     opt_W = problem.variables()[0].value.astype(int)
    
#     # print("Optimal value: {}".format(opt_obj))
#     # print("Solution:\n {}".format(opt_W))
    
#     # matching index Q_idx to P_idx
#     idxQ2P = np.where(opt_W==1)[1]
#     # idxP2Q = np.where(opt_W.T==1)[1]
#     return opt_obj, idxQ2P#, idxP2Q

def alignQ_wrtP(P,Q,idxQ2P,merge=True):
    # K1, K2 = P.shape[1], Q.shape[1]
    
    if merge:        
        aligned_Q = np.zeros_like(P)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        aligned = set()
        extras = list()
        for q_idx in range(Q.shape[1]):
            if idxQ2P[q_idx] in aligned: 
                extras.append(q_idx)
            else:
                aligned.add(idxQ2P[q_idx])
                aligned_Q[:,idxQ2P[q_idx]] = Q[:,q_idx]
        for ie, e in enumerate(extras):
            aligned_Q[:,P.shape[1]+ie] = Q[:,e]
            
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
                # f.write("{}\n".format(" ".join([str(id) for id in idxP2Q])))   

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
         
        align_ILP_res[K] = [[None for _ in ids] for _ in ids]
        cost_ILP_res[K] = [[np.nan for _ in ids] for _ in ids]
        
    for l in range(len(alignments)//2):
        l1 = alignments[2*l].split()
        l2 = alignments[2*l+1].split()
        # l3 = alignments[3*l+2].split()
        
        # cost
        K = int(l1[0])
        ids = k2ids[K]
        i,j = idx2idxinK[int(l1[1])], idx2idxinK[int(l1[2])]
        
        # alignment
        l2 = [int(idx) for idx in l2]
        # l3 = [int(idx) for idx in l3]
        # align_ILP_res[K][i][j] = l3 # P2Q
        align_ILP_res[K][i][j] = [l2.index(idx) for idx in range(K)] # P2Q
        align_ILP_res[K][j][i] = l2 # Q2P
        
        if get_cost:  
            # alignment cost
            P = Q_list[int(l1[1])]
            Q = Q_list[int(l1[2])]
            aligned_idxQ2P = l2
            aligned_Q = np.zeros_like(P)
            for q_idx in range(Q.shape[1]):
                aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
            cost_actual = cost_membership(P,aligned_Q,P.shape[0])
            cost_ILP_res[K][i][j] = cost_actual #float(l1[3])
            cost_ILP_res[K][j][i] = cost_actual #float(l1[3])
        
    if get_cost:
        return align_ILP_res, cost_ILP_res
    else:
        return align_ILP_res

def cd_default(G):
    resolution = 1
    partition_map = community_louvain.best_partition(G,resolution=resolution,random_state=6)
    return partition_map

def cd_custom(G):
    partition_map = {i:0 for i in range(G.number_of_nodes())} 
    return partition_map

def detect_modes(cost_ILP_res,K_range,default_cd_flag,cd_modthre, draw_communities=False, save_path=None, cd_func=cd_default):

    modes_allK_list = dict()
    # adj_mat_allK_list = dict()
    # method_name,method_para = "", None
    # if len(method)>0:
    #     method_name = method[0]
    # if len(method)>1:
    #     method_para = method[1]
    # if len(method)>2:
    #     method_modthre = method[2]
            
    for K in K_range:
    
        # normalize cost matrix and obtain adjacency matrix
        cost_mat = np.array(cost_ILP_res[K])
        adj_mat = get_adj_mat(cost_mat)
        # adj_mat_allK_list[K] = adj_mat
        
        # create graph
        G = nx.from_numpy_matrix(adj_mat)
        G.remove_edges_from(nx.selfloop_edges(G))
        pos = nx.spring_layout(G)
        
        
        ########################################
        # community detection to find modes
        ########################################
        if default_cd_flag:
            partition_map = cd_func(G)
            # print("Running default Louvain for mode detection.")
        else:
            partition_map = cd_custom(G)
            # print("Running customized community detection.")
        if cd_modthre is not None:
            mod_res = community_louvain.modularity(partition_map, G)
            if mod_res<cd_modthre: # quality of community detection is low --> no community structure
                print("K={}: low community detection quality, set to single mode.".format(K))
                partition_map = {i:0 for i in range(G.number_of_nodes())} 
        
        
        # if method_name=="louvain":

        #     resolution = method_para if method_para else 1
        #     partition_map = community_louvain.best_partition(G,resolution=resolution,random_state=6)
        #     modthre = method_modthre if method_modthre else 0.
        #     mod_res = community_louvain.modularity(partition_map, G)
        #     if mod_res<modthre: # quality of community detection is low --> no community structure
        #         print("K={}: low community detection quality, set to single mode.".format(K))
        #         partition_map = {i:0 for i in range(G.number_of_nodes())}
        # else: # no method
        #     print("try louvain in func")
        #     partition_map = cd_func(G)
        #     print("K={}: no community detection algorithm provided, set to single mode".format(K))
        #     partition_map = {i:0 for i in range(G.number_of_nodes())}     

        # # compute the best partition
        # if method_name=="louvain":
        #     resolution = method_para if method_para else 1.
        #     modthre = method_modthre  if method_modthre else 0.
            
        #     coms = algorithms.louvain(G, weight='weight', resolution=resolution, randomize=False)
        #     partition = coms.communities
        #     partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}
        
        #     mod_res = evaluation.newman_girvan_modularity(G,coms).score
        #     if mod_res<modthre: # quality of community detection is low --> no community structure
        #         print("K={}: low community detection quality, set to single mode".format(K))
        #         partition_map = {i:0 for i in range(G.number_of_nodes())}
        
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
    
    
        if draw_communities:
            # draw the graph
            fig, ax = plt.subplots()
            nx.draw(G,pos,with_labels=True,node_color='white',edge_color='white')
            nx.draw_networkx_nodes(G, pos, partition_map.keys(), alpha=0.7, #node_size=40,
                                   node_color=list(partition_map.values()))
            nx.draw_networkx_edges(G, pos, alpha=0.2)
            # save 
            fig.savefig(os.path.join(save_path,'K{}_modes_network_{}.png'.format(K,method_name)))
            plt.close(fig)
    
        ########################################
        # separate modes
        ########################################
        
        modes = defaultdict(list)
        for r in partition_map.keys():
            modes[partition_map[r]].append(r)
        
        modes_allK_list[K] = modes
        
    return modes_allK_list


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
        
        cost_mat = np.array(cost_ILP_res[K])
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
                    # check
                    # cost_recomputed = cost_membership(P,aligned_Q,N_ind)
                    # assert(np.sum(np.abs(cost_recomputed-cost_ILP_res[K][other_self_indices[i]][minC_self_idx]))<0.001)
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

            cost_thre = 5*ct#cost_thre_base+10*np.mean([theo_cost_list[ids[r]],theo_cost_list[ids[leaders[min_cost_idx]]]]) #lc_theo_cost[leaders[min_cost_idx]]
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



# def align_ILP_modes_acrossK(consensusQ_modes,K_range,N,save_path,ILP_acrossK_filename,ind2pop=None):
    
#     K_comb = list(combinations_with_replacement(sorted(K_range,reverse=True),r=2))
#     perm_ILP_acrossK_Q2P = dict()
#     # perm_ILP_acrossK_P2Q = dict()
#     cost_acrossK_ILP = dict()
    
#     f = open(os.path.join(save_path,ILP_acrossK_filename),"w")
    
#     best_alignments = dict()
    
#     for K1,K2 in K_comb:
        
#         best_alignment_idx = None
#         best_alignment_obj = np.inf

#         # labels = list(product(["K{}m{}".format(K1,im) for im,m in enumerate(consensusQ_modes[K1])],\
#         #                       ["K{}m{}".format(K2,im) for im,m in enumerate(consensusQ_modes[K2])]))
#         rijs = list(product(range(len(consensusQ_modes[K1])),range(len(consensusQ_modes[K2]))))
        
#         for i,(ri,rj) in enumerate(rijs):
            
#             Q = consensusQ_modes[K1][ri] # Q is the one with more clusters
#             P = consensusQ_modes[K2][rj]
    
#             opt_obj, idxQ2P = align_ILP(P, Q)
#             if best_alignment_obj>opt_obj:
#                 best_alignment_obj = opt_obj
#                 best_alignment_idx = i
#             f.write("{}#{}-{}#{}:{}\n".format(K2,rj,K1,ri,opt_obj))
#             f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
#             # f.write("{}\n".format(" ".join([str(id) for id in idxP2Q])))
    
#             # retrieve alignment
#             # if "{}#{}-{}#{}".format(K2,rj,K1,ri)=="3#0-3#1":
#             #     print(K1,K2,rijs)
#             #     print("GOOD")
#             #     raise KeyboardInterrupt
#             cost_acrossK_ILP["{}#{}-{}#{}".format(K2,rj,K1,ri)] = opt_obj
#             perm_ILP_acrossK_Q2P["{}#{}-{}#{}".format(K2,rj,K1,ri)] = idxQ2P
#             # perm_ILP_acrossK_P2Q[labels[i]] = idxP2Q
        
#         best_alignments[(K1,K2)] = rijs[best_alignment_idx]
#     f.close()
    
#     # if len(ind2pop)>0:
#     #     f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best.txt"),"w")
#     #     for i_K1 in range(len(K_range)-1):
#     #         K1 = K_range[i_K1+1]
#     #         K2 = K_range[i_K1]
#     #         ri,rj = best_alignments[(K1,K2)]
#     #         Q = consensusQ_modes[K1][ri] # Q is the one with more clusters
#     #         P = consensusQ_modes[K2][rj]
            
#     #         Qbar = get_Qbar([Q],ind2pop)
#     #         Qbar = Qbar[0]
#     #         Pbar = get_Qbar([P],ind2pop)
#     #         Pbar = Pbar[0]
#     #         print(K1,K2,ri,rj,Qbar.shape,Pbar.shape)
    
#     #         opt_obj, idxQ2P = align_IQP(Pbar, Qbar)
#     #         print(opt_obj)
#     #         print(idxQ2P)
#     #         f.write("{}#{}-{}#{}:{}\n".format(K2,rj,K1,ri,opt_obj))
#     #         f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
            
#     #     f.close()
    
#     return perm_ILP_acrossK_Q2P, cost_acrossK_ILP


def align_ILP_modes_acrossK(consensusQ_modes,K_range,N,save_path,ILP_acrossK_filename,ind2pop=None):
    
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
    
            opt_obj, idxQ2P = align_ILP(P, Q)
            if best_alignment_obj>opt_obj:
                best_alignment_obj = opt_obj
                best_alignment_idx = i
            f.write("{}#{}-{}#{}:{}\n".format(K2,rj,K1,ri,opt_obj))
            f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
            # f.write("{}\n".format(" ".join([str(id) for id in idxP2Q])))
    
            # retrieve alignment
            # if "{}#{}-{}#{}".format(K2,rj,K1,ri)=="3#0-3#1":
            #     print(K1,K2,rijs)
            #     print("GOOD")
            #     raise KeyboardInterrupt
            cost_acrossK_ILP["{}#{}-{}#{}".format(K2,rj,K1,ri)] = opt_obj
            perm_ILP_acrossK_Q2P["{}#{}-{}#{}".format(K2,rj,K1,ri)] = idxQ2P
            # perm_ILP_acrossK_P2Q[labels[i]] = idxP2Q
        
        best_alignments[(K1,K2)] = rijs[best_alignment_idx]
    f.close()
    
    # best alignments
    best_ILP_acrossK = list()
    f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best."+ILP_acrossK_filename.split(".")[1]),"w")
    for i_K1 in range(len(K_range_sorted)-1):
        K1 = K_range_sorted[i_K1]
        K2 = K_range_sorted[i_K1+1]
        ri,rj = best_alignments[(K1,K2)]
        bali = "{}#{}-{}#{}".format(K2,rj,K1,ri)
        best_ILP_acrossK.append(bali)
        f.write("{}\n".format(bali))
    f.close()
    
    # if len(ind2pop)>0:
    #     f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best.txt"),"w")
    #     for i_K1 in range(len(K_range)-1):
    #         K1 = K_range[i_K1+1]
    #         K2 = K_range[i_K1]
    #         ri,rj = best_alignments[(K1,K2)]
    #         Q = consensusQ_modes[K1][ri] # Q is the one with more clusters
    #         P = consensusQ_modes[K2][rj]
            
    #         Qbar = get_Qbar([Q],ind2pop)
    #         Qbar = Qbar[0]
    #         Pbar = get_Qbar([P],ind2pop)
    #         Pbar = Pbar[0]
    #         print(K1,K2,ri,rj,Qbar.shape,Pbar.shape)
    
    #         opt_obj, idxQ2P = align_IQP(Pbar, Qbar)
    #         print(opt_obj)
    #         print(idxQ2P)
    #         f.write("{}#{}-{}#{}:{}\n".format(K2,rj,K1,ri,opt_obj))
    #         f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
            
    #     f.close()
    
    return perm_ILP_acrossK_Q2P, cost_acrossK_ILP, best_ILP_acrossK


def load_ILP_acrossK(save_path,ILP_acrossK_filename):

    # load alignments
    f = open(os.path.join(save_path,ILP_acrossK_filename),"r")
    alignments_acrossK = f.readlines()
    f.close()
    
    cost_acrossK_ILP_res = dict()
    # align_acrossK_ILP_res_P2Q = dict()
    align_acrossK_ILP_res_Q2P = dict() # Q has more clusters
    
    n_pairs = (len(alignments_acrossK)-1)//2
    
    for l in range(n_pairs):
        l1 = alignments_acrossK[2*l].split(":")
        l2 = alignments_acrossK[2*l+1].split()
        # l3 = alignments_acrossK[3*l+2].split()
        
        # alignment
        l2 = [int(idx) for idx in l2]
        # l3 = [int(idx) for idx in l3]
    
        cost_acrossK_ILP_res[l1[0]] = l1[1]
        # align_acrossK_ILP_res_P2Q[l1[0]] = l3 # P2Q
        align_acrossK_ILP_res_Q2P[l1[0]] = l2 # Q2P
    
    # best alignments
    f = open(os.path.join(save_path,ILP_acrossK_filename.split(".")[0]+"_best."+ILP_acrossK_filename.split(".")[1]),"r")
    best_ILP_acrossK = f.readlines()
    f.close()
        
    return align_acrossK_ILP_res_Q2P, cost_acrossK_ILP_res, best_ILP_acrossK


def report_stats(modes_allK_list,average_stats,K_range,k2ids,save_path,stats=['cost','Hprime']):

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

