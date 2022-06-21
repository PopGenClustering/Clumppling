# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:15:40 2022

@author: xiranliu
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from TracyWidom import TracyWidom


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


n = 10
W = np.zeros((n,n))
W += np.random.uniform(0,0.1,(n,n))
W[0:n//2,0:n//2] += np.random.uniform(0,0.5,(n//2,n//2))
W[n//2-2:n,n//2-2:n] += np.random.uniform(0,0.5,(n//2+2,n//2+2))
W = np.tril(W, k=-1)
W = W+W.T

has_comm_struc = test_comm_struc(W, alpha = 0.05)
print(has_comm_struc)

# to network
G = nx.from_numpy_matrix(np.matrix(W))
layout = nx.spring_layout(G)
nx.draw(G, layout,  width=list(nx.get_edge_attributes(G,'weight').values()))
# nx.draw_networkx_edge_labels(G, pos=layout)
plt.show()




