# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 19:29:46 2022

@author: xiran
"""

#%% Imports
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
   
# from itertools import product,combinations_with_replacement
# from collections import defaultdict

# set working directory
# os.chdir(os.path.dirname(os.path.realpath(__file__)))


# import community as community_louvain
# from cdlib import algorithms, evaluation
# import networkx as nx
# import markov_clustering as mc # MCL for mode detection

# import networkx.algorithms.community as nx_comm
# import markov_clustering as mc # MCL for mode detection

# import community.community_lc_cost_thre

# from scipy.spatial.distance import cdist
# import cvxpy as cp


from clumppling.func_utils import *
from clumppling.func_main import *
from clumppling.func_plotting import *
from clumppling.params import Params
import argparse


def main(args):

    input_path = args.input_path
    output_path = args.output_path
    prj_type = args.prj_type
    
    if not prj_type in ["structure","fastStructure","admixture","generalQ"]:
        sys.exit("ERROR: Input project type isn't supported. \nPlease specify prj_type as one of the following: structure, admixture, fastStructure, and generalQ.")
    
  
    params = Params(input_path,output_path,prj_type)
    
    if args.default_cd:
        params.default_cd = args.default_cd
    if args.cd_mod_thre:
        params.cd_mod_thre = args.cd_mod_thre
    if args.reorder_inds:
        params.reorder_inds = args.reorder_inds
    if args.custom_cmap:
        params.custom_cmap = args.custom_cmap
        params.cmap = args.cmap.split()
    
    if os.path.exists(params.input_path+".zip"):
        shutil.unpack_archive(params.input_path+".zip",params.input_path)
    if not os.path.exists(params.input_path):
        sys.exit("Input file doesn't exist.")
    
    #%% Set-up 
    tot_tic = time.time()
    
    # cmap = cm.get_cmap('Spectral') # colormap for plotting clusters
    # cmap_modes = cm.get_cmap('tab10') # colormap for plotting mode network
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    if params.lc_flag:
        md_method = "LC_{}".format("adaptive_{}".format(params.lc_cost_thre) if params.adaptive_thre_flag else params.lc_cost_thre) 
    else:
        md_method = "{}_{}".format("default" if params.default_cd else "custom",params.cd_mod_thre)
    
    if params.Qbar_flag:
        save_path = os.path.join(params.output_path, "Qbar_"+md_method) 
    else:
        save_path = os.path.join(params.output_path, md_method)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    output_f = os.path.join(save_path,'output.txt')
    handlers = [logging.FileHandler(output_f, 'w'), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format='', handlers = handlers)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True

    recode_path = os.path.join(params.output_path,"data")
    if not os.path.exists(recode_path) or len(os.listdir(recode_path))==0:
        # remove empty directory if exist
        if os.path.exists(recode_path):
            os.rmdir(recode_path)
        
        os.makedirs(recode_path) 
        tic = time.time()
        
        if params.prj_type =="structure":
            recode_struture_files(params.input_path,recode_path)
        elif params.prj_type =="fastStructure":
            recode_faststruture_files(params.input_path,recode_path)
        elif params.prj_type =="admixture" or params.prj_type =="generalQ":
            recode_admixture_files(params.input_path,recode_path)
        else:
            sys.exit("Input project type isn't supported. Please specify one of the following: structure, admixture, fastStructure, and generalQ.")
        
        toc = time.time()
        logging.info("time to recode structure files: %s",toc-tic)
    
    
    #%% Loading membership data
    tic = time.time()
    
    N, R, Q_list, K_list = load_Q(recode_path,reorder_inds=params.reorder_inds)
    if prj_type == "structure":
        ind2pop, pop_n_ind = load_ind(recode_path)
    else:
        ind2pop = None
    
    K_list = np.array(K_list)
    K_range = np.sort(np.unique(K_list))
    max_K = max(K_range)
    k2ids = {k:np.where(K_list==k)[0] for k in K_range}
    idx2idxinK = [i-np.where(K_list==K)[0][0] for i,K in enumerate(K_list)]
    
    # colormap
    if params.custom_cmap:
        if len(params.cmap) < max_K:
            logging.info("The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            params.cmap.extend(params.cmap)
        cmap = cm.colors.ListedColormap(params.cmap)
        logging.info("using custom cmap ...")
    else:
        np.random.seed(9999)
        cmap = cm.get_cmap('nipy_spectral') # colormap for plotting clusters
        cmap = cmap(np.linspace(0, 1, max_K))
        np.random.shuffle(cmap)
        cmap = cm.colors.ListedColormap(cmap)
    
    # color of mode network layers
    cmap_modes = cm.get_cmap('tab10') # colormap for plotting mode network
    layer_color = {K:cmap_modes(i) for i,K in enumerate(K_range)}
    
    toc = time.time()
    logging.info("time to load and set up replicates: %s",toc-tic)
    
    #%% Alignment within-K and mode detection
    
    tic = time.time()
    ILP_modes_filename = "ILPmodes.txt"
    
    if params.lc_flag:
        ## Leader Clustering
        if params.adaptive_thre_flag:
            modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats = align_leader_clustering_adaptive(lc_cost_thre,Q_list,K_range,N,k2ids,ind2pop, pop_n_ind,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
        else:
            modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats = align_leader_clustering(lc_cost_thre,Q_list,K_range,N,k2ids,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
        
    else:
        ## ILP over all pairs within K 
        # write to file then load from file
        ILP_withinK_filename = "ILPaligned.txt"
        
        if params.Qbar_flag:
            ind2pop, pop_n_ind = load_ind(recode_path)
            Qbar_list = get_Qbar(Q_list,ind2pop)
            if not os.path.exists(os.path.join(params.output_path,ILP_withinK_filename)):
                # write
                align_ILP_withinK(ILP_withinK_filename,params.output_path,Qbar_list,K_range,k2ids)
            else:
                logging.info("pairwise alignment (Qbar) already exist")
        else:
            if not os.path.exists(os.path.join(params.output_path,ILP_withinK_filename)):
                # write
                align_ILP_withinK(ILP_withinK_filename,params.output_path,Q_list,K_range,k2ids)
            else:
                logging.info("pairwise alignment already exist")
        # load
        align_ILP_res, cost_ILP_res = load_ILP_withinK(Q_list,ILP_withinK_filename,output_path,K_range,k2ids,idx2idxinK,get_cost=True)
        
        ## Mode detection 
        modes_allK_list = detect_modes(cost_ILP_res,K_range,params.default_cd,params.cd_mod_thre,draw_communities=params.plot_flag_community_detection, save_path=save_path)
           
        # extract representative/consensus replicates
        rep_modes,repQ_modes,meanQ_modes,average_stats = extract_modes(Q_list,modes_allK_list,align_ILP_res,cost_ILP_res,K_range,N,k2ids,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
    
    ## Within-K performance
    # stats of alignments: cost and hprime (similarity)
    weighted_stats = report_stats(modes_allK_list,average_stats,K_range,k2ids,save_path=save_path,stats=['cost','Hprime'])
    
    
    toc = time.time()
    logging.info("time to align replicates and detect modes: %s", toc-tic)
    logging.info("(method used: {})".format(md_method))
    
    #%% Alignment across-K
        
    # use representative membership as the consensus
    ILP_acrossK_filename = "ILP_acrossK_repQ.txt"
    
    if not os.path.exists(os.path.join(save_path,ILP_acrossK_filename)):
        # write
        tic = time.time()
        repQ_acrossK_Q2P, repQ_acrossK_cost, repQ_best_ILP_acrossK = align_ILP_modes_acrossK(repQ_modes,K_range,N,save_path,ILP_acrossK_filename,ind2pop=ind2pop)
        toc = time.time()
        logging.info("time to ILP over all modes (repQ) across-K: %s", toc-tic)
    else:
        # load alignments 
        repQ_acrossK_Q2P, repQ_acrossK_cost, repQ_best_ILP_acrossK = load_ILP_acrossK(save_path,ILP_acrossK_filename)
        logging.info("across-K alignment (repQ) already exist")
    

    # use average membership as the consensus
    ILP_acrossK_filename = "ILP_acrossK_meanQ.txt"
    
    if not os.path.exists(os.path.join(save_path+ILP_acrossK_filename)):
        # write
        tic = time.time()
        meanQ_acrossK_Q2P, meanQ_acrossK_cost, meanQ_best_ILP_acrossK = align_ILP_modes_acrossK(meanQ_modes,K_range,N,save_path,ILP_acrossK_filename)
        toc = time.time()
        logging.info("time to ILP over all modes (meanQ) across-K: %s", toc-tic)
    else:
        meanQ_acrossK_Q2P, meanQ_acrossK_cost, meanQ_best_ILP_acrossK = load_ILP_acrossK(save_path,ILP_acrossK_filename)
        logging.info("across-K alignment (meanQ) already exist")
    
    #%% Visualization of alignment within-K
    tic = time.time()
    ##  Alignment of each mode
    if params.plot_flag_all_within_mode:
        for K in K_range: #[1:2]:
            modes = modes_allK_list[K]
            for m in modes:
                plot_name = "K{}_mode{}_meanQ.png".format(K,m)
                plot_aligned(K,m,Q_list,modes,align_ILP_res,rep_modes,meanQ_modes[K][m],max_K,k2ids,idx2idxinK,save_path=save_path,plot_name=plot_name,cmap=cmap)
    
    ## Alignment between modes within-K
    if params.plot_flag_mode_within_K:
        plot_file_name_suffix = "modes"
        for K in K_range:
            plot_withinK_modes(K,max_K,meanQ_modes,meanQ_acrossK_Q2P,save_path,plot_file_name_suffix,cmap=cmap)
            
    toc = time.time()
    logging.info("time to plot alignment within-K: %s", toc-tic)
    
    #%% Visualization of alignment across-K
    tic = time.time()
    ## Average membership multipartite graph
    if params.plot_flag_mode_across_K_multipartite:
        plot_file_name = "acrossK_meanQ_cost.png"
        title = "Mode (average membership) across K"
        G = plot_acrossK_multipartite(K_range,modes_allK_list,meanQ_modes,meanQ_acrossK_cost,layer_color,title,save_path,plot_file_name)
        
        plot_file_name = "acrossK_repQ_cost.png"
        title = "Mode (representative membership) across K"
        G = plot_acrossK_multipartite(K_range,modes_allK_list,repQ_modes,repQ_acrossK_cost,layer_color,title,save_path,plot_file_name)
    
    ## Alignment of all modes
    if params.plot_flag_all_modes:
        plot_name = "all_modes.png" 
        plot_all_modes(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_best_ILP_acrossK,save_path,plot_name,cmap=cmap)    
    
    ## Optimal alignment across-K in chains
    if params.plot_flag_mode_across_K_chains:
        plot_file_name_suffix = "meanQ_alignment_chain_merged"
        best_alignment_chains = plot_acrossK_chains(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_acrossK_cost,save_path,plot_file_name_suffix,cmap=cmap,merge_cluster=True)
        
        plot_file_name_suffix = "meanQ_alignment_chain"
        best_alignment_chains = plot_acrossK_chains(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_acrossK_cost,save_path,plot_file_name_suffix,cmap=cmap,merge_cluster=False)
    
    toc = time.time()
    logging.info("time to plot alignment across-K: %s", toc-tic)
    
    
    tot_toc = time.time()
    logging.info("Total Time: %s", tot_toc-tot_tic)
    
    # zip all files
    shutil.make_archive(params.output_path, "zip", params.output_path)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str, required=True)
    parser.add_argument('--output_path', type=str, required=True)
    parser.add_argument('--prj_type', type=str, required=True)
    
    optional_arguments = [['default_cd','bool','whether to use default community detection method (Louvain) to detect modes'],
                          ['cd_mod_thre','float','the modularity threshold for community detection'],
                          ['reorder_inds','bool','whether to reorder individuals based on memberships from the last Q with largest K'],
                           ['custom_cmap','bool','whether to use customized colormap'],
                           ['cmap','str','user-specified colormap as a list of colors (in hex code) in a space-delimited string']]
    
    for opt_arg in optional_arguments: 
        parser.add_argument('--{}'.format(opt_arg[0]), type=getattr(builtins, opt_arg[1]), required=False, help=opt_arg[2])

    
    # parser.add_argument('--default_cd', type=bool, required=False)
    # parser.add_argument('--cd_mod_thre', type=float, required=False)
    # parser.add_argument('--reorder_inds', type=bool, required=False)
    args = parser.parse_args()
    main(args)