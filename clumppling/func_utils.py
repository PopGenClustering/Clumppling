# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:25:05 2022

@author: xiranliu
"""

import numpy as np
import pandas as pd
import os
import sys
# import matplotlib.pyplot as plt
import time
from itertools import product,combinations_with_replacement
from collections import defaultdict

from scipy import special

# import community as community_louvain
# from cdlib import algorithms, evaluation
# import community as community_louvain
# import networkx as nx
# import markov_clustering as mc # MCL for mode detection

# import networkx.algorithms.community as nx_comm
# import markov_clustering as mc # MCL for mode detection

# import community.community_louvain
# import matplotlib.cm as cm


# import scipy.stats as ss

# from scipy.spatial.distance import cdist
# import cvxpy as cp

#%% Helper Functions

def cost_membership(P,Q,N_ind):
    return np.linalg.norm(P-Q, ord="fro")**2/(2*N_ind)

def clumppGprime(P,Q,N_ind):
#     return 1-np.linalg.norm(P-Q, ord="fro")/np.sqrt(2*N_ind)
    return 1-np.sqrt(cost_membership(P,Q,N_ind))

def C2Gprime(C):
    return 1-np.sqrt(C)

def avg_tril(A):
    return np.sum(np.tril(A, k=-1))/(A.shape[0]*(A.shape[0]-1)/2)


#%% Data Prep Functions

# def reorder_Q_inds(Q_ind):
#     Q_ind_reorder = list(np.hstack((Q_ind,np.expand_dims(np.arange(Q_ind.shape[0],dtype='object'),axis=-1))))
    
#     Q_ind_reorder.sort(key=lambda x:(x[0], x[1] ),reverse=True)
#     Q_ind_reorder = np.array(Q_ind_reorder)
#     reorder_idx = Q_ind_reorder[:,-1].astype('int32')
#     return reorder_idx


def recode_files(data_path,recode_path,input_format):


    files = os.listdir(data_path)
    indivq_ext = False
    if input_format == "structure":
        Q_files = [i for i in files if i.endswith('_f')]
    elif input_format == "fastStructure":
        Q_files = [i for i in files if i.endswith('.meanQ')]
    elif input_format == "admixture" or input_format == "generalQ":
        Q_files = [i for i in files if i.endswith('.Q')]
        if len(Q_files)==0: # indivq files
            Q_files = [i for i in files if i.endswith('.indivq')]
            indivq_ext = True
    else:
        sys.exit("ERROR: Please specify one of the following for input_format: structure, admixture, fastStructure, and generalQ.")

    if len(Q_files)==0:
        # check if no input data of the right form 
        sys.exit("ERROR: No input files detected. Please double check input_path and input_format.")

 
    R = len(Q_files)
    
    if input_format =="structure":
        for r in range(R):
            Q_file = Q_files[r]
            with open(os.path.join(data_path,Q_file)) as file:
                lines = file.readlines()
            
            # get only the membership coefficients part
            res_start = lines.index("Inferred ancestry of individuals:\n")
            lines = lines[(res_start+1):]
            res_end = lines.index("\n")
            lines = lines[:res_end]
            
            # extract membership matrix
            Q = [[float(i) for i in l.split()[5:]] for l in lines[1:]]
            Q = np.array(Q)
            
            # write to file
            np.savetxt(os.path.join(recode_path,'rep_{}.Q'.format(r)), Q, delimiter=' ')

    elif indivq_ext:
        for r in range(R):
            Q_file = Q_files[r]
            Q = np.loadtxt(os.path.join(data_path,Q_file),dtype=str)[:,5:].astype(float)
            np.savetxt(os.path.join(recode_path,'rep_{}.Q'.format(r)), Q, delimiter=' ')

    else: # if input_format =="admixture" or "fastStructure" or "generalQ"
        for r in range(R):
            Q_file = Q_files[r]
            Q = np.loadtxt(os.path.join(data_path,Q_file)).astype(float)

            # with open(os.path.join(data_path,Q_file)) as file:
            #     lines = file.readlines()
                
            # # extract membership matrix
            # Q = [[float(i) for i in l.split()] for l in lines]
            # Q = np.array(Q)
            
            # write to file
            np.savetxt(os.path.join(recode_path,'rep_{}.Q'.format(r)), Q, delimiter=' ')

    if input_format =="structure":
        # extract individual info
        df_ind_info = pd.DataFrame([l.split()[0:4] for l in lines[1:]])
        df_ind_info.set_index([0],inplace=True)
        header = lines[0].split()[:3]
        header[2] = header[2].split(":")[0]
        df_ind_info.columns = header
        df_ind_info.to_csv(os.path.join(recode_path,'ind_info.txt'), sep=' ')

    return Q_files

    
def load_Q(recode_path,file_list=None):

    Q_list = list()
    K_list = list()
    # pop_n_ind = None
    
    Q_files = [i for i in os.listdir(recode_path) if i.endswith('.Q')]
       
    R = len(Q_files)
    if R==0:
        # sanity check if data directory is empty
        sys.exit("ERROR: no Q files detected. Please double check input_path and input_format.")
    

    for r in range(R):
        if not file_list:
            Q = np.loadtxt(os.path.join(recode_path,Q_files[r]))
        else:
            Q = np.loadtxt(os.path.join(recode_path,'rep_{}.Q'.format(r)))
        
        K = Q.shape[1]
        if r==0:
            N = Q.shape[0]
        else:
            assert(Q.shape[0]==N)
            
        Q_list.append(Q)
        K_list.append(K)
    
    if not file_list:
        file_list = Q_files
        
    # reorder by consecutive K values   
    sored_idx = list(np.argsort(K_list))
    K_list = [K_list[i] for i in sored_idx]
    Q_list = [Q_list[i] for i in sored_idx]
    num_zeros_in_all_Q = np.sum([np.sum(Q==0) for Q in Q_list])
    file_list = [file_list[i] for i in sored_idx]
                   
    N = Q.shape[0]

    

    # old file to Q idx
    with open(os.path.join(recode_path,'idx2file.txt'),'w') as file:
        for r,f_name in enumerate(file_list):
            file.write('{}\t{}\n'.format(r,f_name))
    
    return N, R, Q_list, K_list, file_list


def load_ind(recode_path):
        
    info = pd.read_csv(os.path.join(recode_path,'ind_info.txt'),delimiter=" ",index_col=0)
    # get populization sizes
    popIDs = info["Pop"].unique()
    pop_n_ind = {}
    for p in popIDs:
        pop_n_ind[p] = sum(np.array(info["Pop"]==p))
    ind2pop = info["Pop"].values
    
    return ind2pop, pop_n_ind

# def get_Qbar(Q_list,ind2pop):
#     ind_pop_first = [i+1 for i,ind in enumerate(ind2pop[:-1]) if ind!=ind2pop[i+1] ]    
#     ind_pop_first = [0] + ind_pop_first + [len(ind2pop)]

#     Qbar_list = list()
#     for Q in Q_list:
#         Q_pop = [np.sum(Q[ind_pop_first[i]:ind_pop_first[i+1],:],axis=0,keepdims=True)/(ind_pop_first[i+1]-ind_pop_first[i])*(ind_pop_first[i+1]-ind_pop_first[i]) for i in range(len(ind_pop_first)-1)]
#         Q_pop = np.array(Q_pop).squeeze()
#         Qbar_list.append(Q_pop)
        
#     return Qbar_list

#%% Dirichlet Functions
def digamma_inv(y,n_iter=5):
    # Inverse of the digamma function (real positive arguments only).
    #
    gamma = -special.psi(1)
    t_thre = -2.22
    if y > t_thre:
        x = np.exp(y) + 0.5
    else:
        x = 1.0 / (-y - gamma)

    for i_iter in range(n_iter):
        x = x-(special.psi(x)-y)/special.polygamma(1, x)

    return x

def initial_guess(Q):
    Eq1 = np.mean(Q[:,0])
    Eq1sqr = np.mean(Q[:,0]**2)
    frac = Q.sum(axis=0)/Q.sum()
    if Q.shape[0]==1:
        return frac
    denom = (Eq1sqr-Eq1**2)
    if (Eq1sqr-Eq1**2)!=0:
        a = np.multiply(frac,(Eq1-Eq1sqr)/(Eq1sqr-Eq1**2))
    else:
        a = frac
    if np.all(np.abs(a)<1e-9) or np.all(np.abs(a)>1e9):
        a = frac
    return a
  
def fixed_point(Q, a0, n_iter = 10):
    K = len(a0)
    N = Q.shape[0]
    if N==1:
        return a0
    logq_bar = np.sum(np.log(Q),axis=0)/N
    a = a0
    a_next = a
    for i_iter in range(n_iter):
        a_sum = np.sum(a)
        for k in range(K):
            a_next[k] = digamma_inv(special.psi(a_sum)+logq_bar[k])
        if np.sum(np.abs(a_next-a))<1e-3:
            a = a_next
            break
        a = a_next
#     print("#iterations:{}".format(i_iter+1))
    return a

def repdist0(a):
    a0 = np.sum(a)
    return 4*np.sum(np.tril(np.outer(a,a),-1))/((a0+1)*(a0**2))

def get_theoretical_cost(Q,pop_n_ind,ind2pop):
    a_vec_list = []
    theo_cost = []
    theo_cost_ws = 0
    for p in pop_n_ind.keys():
    
        Q_pop = Q[ind2pop==p,:]
        a0 = initial_guess(Q_pop)
        a = fixed_point(Q_pop, a0) # estimated parameter
        if np.isnan(a[0]):
            break
        a_vec_list.append(a)
        theo_cost.append(repdist0(a))
        theo_cost_ws += pop_n_ind[p]*repdist0(a)
        
    theo_cost_ws *= 2/len(ind2pop)

    return theo_cost_ws


def get_Dir_params(Q,pop_n_ind,ind2pop):

    # infer Dir parameters
    a_vec_list = []
    
    for p in pop_n_ind.keys():
    
        Q_pop = Q[ind2pop==p,:]
        if np.sum(Q_pop==0)>0:
            Q_pop += 1e-6
        a0 = initial_guess(Q_pop)
        a = fixed_point(Q_pop, a0) # estimated parameter
        if np.isnan(a[0]):
            break
        a_vec_list.append(a)

    return a_vec_list

def cost_batch(a,perms,coords):
    
    all_entry = (a[:,None]-a[None,:])*a[:, np.newaxis]
    
    flat_idx = np.concatenate((np.asarray(all_entry.shape[1:])[::-1].cumprod()[::-1],[1])).dot(coords)
    a0 = np.sum(a)
    c = np.take(all_entry,flat_idx).reshape(perms.shape)
    c = np.sum(c,axis=1)/a0**2
    return c

def compute_cost_total(a_list,n_ind_list,perms):
    N = sum(n_ind_list)
    K = len(a_list[0])
    coords = np.vstack([np.tile(np.arange(K), len(perms)),perms.flatten()])
    the_C_pop = []
    for a in a_list:
        c = cost_batch(a,perms,coords)
        the_C_pop.append(c)
    the_C_pop = np.array(the_C_pop)
    the_C_total = np.sum(the_C_pop*np.array(n_ind_list)[:, np.newaxis]/N,axis=0)
    return the_C_total


def sort_ignore_tie(l):
    l = np.array(l)
    sorted_idx = l.argsort()
    sorted_val = l[sorted_idx]
    return sorted_idx, sorted_val



#%% Mode Detection Helpers
def get_adj_mat(cost_mat):
    adj_mat = 1-cost_mat/np.nanmax(cost_mat)
    adj_mat = np.nan_to_num(adj_mat,copy=True,nan=0)
    adj_mat = adj_mat + np.diag(np.ones(adj_mat.shape[0]))
    return adj_mat


def write_modes_to_file(file_name,K_range,N_ind,modes_allK_list,meanQ_modes):
    f = open(file_name,"w")
    for i_K, K in enumerate(K_range):
        
        f.write("K={}|{}|{}\n".format(K,i_K,N_ind))
        for mode_idx in range(len(meanQ_modes[K])):
            indices_in_mode = " ".join([str(item) for item in modes_allK_list[K][mode_idx]])
            meanQ = meanQ_modes[K][mode_idx]
            meanQ_str = ' '.join(str(np.round(item,4)) for innerlist in meanQ for item in innerlist)
            f.write("mode|{}|{}\n".format(mode_idx, indices_in_mode))
            f.write("{}\n".format(meanQ_str))
        f.write("\n")
    
    f.close()


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
