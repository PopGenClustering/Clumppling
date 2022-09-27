# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:45:25 2022

@author: xiranliu
"""

import numpy as np
# import pandas as pd
import os
import sys
# import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import logging
import shutil
import builtins
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

from clumppling.func_utils import load_Q, load_ind
from clumppling.func_main import align_ILP
from clumppling.func_plotting import plot_membership
import argparse


def load_Q_multiple(input_base_path,input_names):
    # load mode Q   
    N_all = list()
    R_all = list()
    Q_all = list()
    K_all = list()

    for input_name in input_names:
        Q_path = os.path.join(input_base_path, input_name,"modes_Q")
        N, R, Q_list, K_list, file_list = load_Q(Q_path,ignore_recode_name=True,file_list=None)
        N_all.append(N)
        R_all.append(R)
        Q_all.append(Q_list)
        K_all.append(K_list)
        
    return N_all, R_all, Q_all, K_all
    

def align_multiple(input_names, Q_all, K_all, output_path, cmap, plot_separate=False):

        
    K_range = np.sort(list(set().union(*[set(K) for K in K_all])))
    max_K = np.max(K_range)
    
    num_input = len(input_names)
    
    # align modes 
    cost_ILP_res = dict()
    align_ILP_res = dict()
    
    f = open(os.path.join(output_path,"ILP_alignment.txt"),"w")
    
    for K in K_range[::-1]:
        Qrep_list = list()
        # align each K (to input 0)
        for q in range(num_input):
            Q = Q_all[q][np.where(np.array(K_all[q])==K)[0][0]]
            Qrep_list.append(Q)
        
        cost_ILP_res[K] = [[np.nan for _ in input_names] for _ in input_names]  
        align_ILP_res[K] = [[None for _ in input_names] for _ in input_names]    
        
        for i_p in range(num_input):
            P = Qrep_list[i_p]

            for i_q in range(num_input):
                Q = Qrep_list[i_q]
                            
                opt_obj, idxQ2P = align_ILP(P, Q)
                
                f.write("{} {} {} {}\n".format(K,i_p,i_q,opt_obj))
                f.write("{}\n".format(" ".join([str(id) for id in idxQ2P])))
                
                cost_ILP_res[K][i_p][i_q] = opt_obj
                align_ILP_res[K][i_p][i_q] = idxQ2P
    
    f.close()
    
    # plot alignment
    
    K_maxcnt = defaultdict(lambda: float('-inf'))
    for d in [Counter(K) for K in K_all]:
        for k, v in d.items():
            K_maxcnt[k] = max(K_maxcnt[k], v)
    K_cnt_list = [K_maxcnt[k] for k in K_range]
    tot_K_cnt = np.sum(K_cnt_list)
    sep_K_cnt = np.cumsum(K_cnt_list)
    sep_K_cnt = np.insert(sep_K_cnt, 0, 0)
    
    if plot_separate:

        i_p = 0
        num_K = len(K_all[i_p])
        fig, axes = plt.subplots(num_K,1,figsize=(20,2*(num_K)))
        
        fig_idx = 0
        for i_K,K in enumerate(K_range):
            for idx in np.where(np.array(K_all[i_p])==K)[0]:
                P = Q_all[i_p][idx]
                ax = axes[fig_idx]
                plot_membership(ax,P,max_K,cmap,"K={}".format(K))
                fig_idx += 1 
        fig.savefig(os.path.join(output_path,"aligned_{}.pdf".format(input_names[i_p])), bbox_inches='tight', dpi=20)
        plt.close(fig)
        
        for i_q in range(1,num_input):
            num_K = len(K_all[i_q])
            fig, axes = plt.subplots(num_K,1,figsize=(20,2*(num_K)))
            
            fig_idx = 0
            for i_K,K in enumerate(K_range):
                for idx in np.where(np.array(K_all[i_q])==K)[0]:
                    Q = Q_all[i_q][idx]
                    idxQ2P = align_ILP_res[K][i_p][i_q]
                    aligned_Q = np.zeros_like(Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
                    ax = axes[fig_idx]
                    plot_membership(ax,aligned_Q,max_K,cmap,"K={}".format(K))
                    fig_idx += 1
            fig.savefig(os.path.join(output_path,"aligned_{}.pdf".format(input_names[i_q])), bbox_inches='tight', dpi=20)
            plt.close(fig)
                
    else:
    
        fig, axes = plt.subplots(tot_K_cnt,num_input,figsize=(30,2*(tot_K_cnt)))
        
        for i_K,K in enumerate(K_range):
            i_p = 0
            k_idx = 0
            for idx in np.where(np.array(K_all[i_p])==K)[0]:
                P = Q_all[i_p][idx]
                ax = axes[k_idx+sep_K_cnt[i_K],i_p]
                plot_membership(ax,P,max_K,cmap,"K={}".format(K) if k_idx==0 else "")
                k_idx += 1
            while k_idx<K_maxcnt[K]:
                ax = axes[k_idx+sep_K_cnt[i_K],i_p]
                ax.axis('off')
                k_idx += 1
            
            for i_q in range(1,num_input):
                k_idx = 0
                for idx in np.where(np.array(K_all[i_q])==K)[0]:
                    Q = Q_all[i_q][idx]
                    idxQ2P = align_ILP_res[K][i_p][i_q]
                    aligned_Q = np.zeros_like(Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
                    ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                    plot_membership(ax,aligned_Q,max_K,cmap,"")
                    k_idx += 1
                while k_idx<K_maxcnt[K]:
                    ax = axes[k_idx+sep_K_cnt[i_K],i_q]
                    ax.axis('off')
                    k_idx += 1
            
        for i_p in range(num_input):
            ax = axes[0,i_p]
            ax.set_title(input_names[i_p], fontsize=18)
            
        fig.savefig(os.path.join(output_path,"aligned_all.pdf"), bbox_inches='tight', dpi=20)
        plt.close(fig)


def main(args):
    
    input_base_path = args.input_base_path
    output_path = args.output_path

    if args.plot_separate:
        if args.plot_separate == "Y":
            plot_separate = True
        else: 
            plot_separate = False
    else:
        plot_separate = False

    if args.custom_cmap:
        custom_cmap = True if args.custom_cmap=="Y" else False    
    if args.cmap:
        cmap = args.cmap.split()
    else:
        cmap = []
    
    # create output directory
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    if os.path.exists(output_path+".zip"):
        os.remove(output_path+".zip")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    tot_tic = time.time()
    
    output_f = os.path.join(output_path,'output.log')
    handlers = [logging.FileHandler(output_f, 'w'), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format='', handlers = handlers)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True
    
    logging.info("========== Starting Clumppling (different models) ... ==========")
    
    # determine input directories
    input_names = [f for f in os.listdir(input_base_path) if os.path.isdir(os.path.join(input_base_path, f))]
    logging.info("Input Q files: {}".format(", ".join(input_names)))
    logging.info("Plotting figures separately: {}".format("Yes" if plot_separate else "No"))

    # load files
    tic = time.time()
    logging.info("---------- Loading files ...")
    N_all, R_all, Q_all, K_all = load_Q_multiple(input_base_path,input_names)
    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)
    
    # alignment
    tic = time.time()
    logging.info("---------- Aligning and plotting ...")
    
    K_range = np.sort(list(set().union(*[set(K) for K in K_all])))
    max_K = np.max(K_range)

    # set colormap
    if custom_cmap and len(cmap)>0:
        while len(cmap) < max_K:
            logging.info(">>>The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            cmap.extend(cmap)
        cmap = cm.colors.ListedColormap(cmap)
    else:
        if custom_cmap:
            logging.info(">>>Custom colormap is not provided. Use the deafult colormap.")
        np.random.seed(999)
        cmap = cm.get_cmap('Spectral') # colormap for plotting clusters
        cmap = cmap(np.linspace(0, 1, max_K))
        np.random.shuffle(cmap)
        cmap = cm.colors.ListedColormap(cmap)


    # align and plot
    align_multiple(input_names, Q_all, K_all, output_path, cmap, plot_separate=plot_separate)
    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)    
    
    # zip all files
    logging.info("---------- Zipping files ...")
    shutil.make_archive(output_path, "zip", output_path)

    tot_toc = time.time()
    logging.info("======== Total Time: %.3fs ========", tot_toc-tot_tic)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input_base_path', type=str, required=True)
    required.add_argument('-o', '--output_path', type=str, required=True)
    optional.add_argument('--plot_separate', type=str, required=False, help='Y/N: whether plot alignment separately for each model')
    optional.add_argument('--custom_cmap', type=str, required=False, help='Y/N: whether to use customized colormap')
    optional.add_argument('--cmap', type=str, required=False, help='user-specified colormap as a list of colors (in hex code) in a space-delimited string')

    args = parser.parse_args()
    main(args)