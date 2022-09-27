# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 19:29:46 2022

@author: xiran
"""

#%% Imports
import numpy as np
import os
import sys
import matplotlib.cm as cm
import time
import logging
import shutil
import builtins   


from clumppling.func_utils import *
from clumppling.func_main import *
from clumppling.func_plotting import *
from clumppling.params import Params
import argparse

# import xml.dom.minidom


def main(args):

    input_path = args.input_path
    output_path = args.output_path
    params_path = args.params_path
    input_format = args.input_format
    
    # sanity check for arguments
    if os.path.exists(input_path+".zip"):
        shutil.unpack_archive(input_path+".zip",input_path)
    if not os.path.exists(input_path):
        sys.exit("ERROR: Input file doesn't exist.")
    if not os.path.exists(params_path):
        sys.exit("ERROR: Parameter file doesn't exist.")
    if not input_format in ["structure","fastStructure","admixture","generalQ"]:
        sys.exit("ERROR: Input data format isn't supported. \nPlease specify input_format as one of the following: structure, admixture, fastStructure, and generalQ.")
    
    params = Params(input_path,output_path,params_path,input_format)

    if args.default_cd:
        params.default_cd = True if args.default_cd=="Y" else False        
    if args.cd_mod_thre:
        params.cd_mod_thre = args.cd_mod_thre if args.cd_mod_thre!=-1 else -1.0
 
    if args.custom_cmap:
        params.custom_cmap = True if args.custom_cmap=="Y" else False  
    if args.cmap:
        params.cmap = args.cmap.split()
    else:
        params.cmap = []
    
    
    
    #%% Set-up 
    tot_tic = time.time()
    
    # cmap = cm.get_cmap('Spectral') # colormap for plotting clusters
    # cmap_modes = cm.get_cmap('tab10') # colormap for plotting mode network
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    if os.path.exists(output_path+".zip"):
        os.remove(output_path+".zip")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    save_path = output_path
        
    # if not os.path.exists(os.path.join(save_path,"modes_Q")):
    os.makedirs(os.path.join(save_path,"modes_Q"))
    
    output_f = os.path.join(save_path,'output.log')
    handlers = [logging.FileHandler(output_f, 'w'), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format='', handlers = handlers)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True

    disp = params.display()
    logging.info(disp)

    logging.info("========== Starting Clumppling ... ==========")
    # logging.info("Mode detection method: {}".format("community detection" if not params.lc_flag else "leader clustering"))
    # if params.lc_flag:
    #     logging.info("---Mode detection method: LC \n---LC cost threshold: {}".format("adaptive_{}".format(params.lc_cost_thre) if params.adaptive_thre_flag else params.lc_cost_thre))
    # else:
    #     logging.info("---Mode detection method: {} \n---Modularity threshold for mode detection: {}".format("default" if params.default_cd else "custom",params.cd_mod_thre))
    # logging.info("==========")

    recode_path = os.path.join(params.output_path,"data")
    
    os.makedirs(recode_path) 
    tic = time.time()
    logging.info("---------- Processing input data files and checking arguments ...")
    file_list = recode_files(params.input_path,recode_path,params.input_format)
    
    
    #%% Loading membership data    
    N, R, Q_list, K_list, file_list = load_Q(recode_path,file_list=file_list)
   
    if input_format == "structure":
        ind2pop, pop_n_ind = load_ind(recode_path)
        logging.info(">>>Extract individual information from STRUCTURE files.")
    else:
        ind2pop = None
    
    K_list = np.array(K_list)
    K_range = np.sort(np.unique(K_list))
    max_K = max(K_range)
    k2ids = {k:np.where(K_list==k)[0] for k in K_range}
    idx2idxinK = [i-np.where(K_list==K)[0][0] for i,K in enumerate(K_list)]
    
    # set colormap
    if params.custom_cmap and len(params.cmap)>0:
        while len(params.cmap) < max_K:
            logging.info(">>>The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            params.cmap.extend(params.cmap)
        cmap = cm.colors.ListedColormap(params.cmap)
    else:
        if params.custom_cmap:
            logging.info(">>>Custom colormap is not provided. Use the deafult colormap.")
        np.random.seed(999)
        cmap = cm.get_cmap('Spectral') # colormap for plotting clusters
        cmap = cmap(np.linspace(0, 1, max_K))
        np.random.shuffle(cmap)
        cmap = cm.colors.ListedColormap(cmap)
    
    # color of mode network layers
    cmap_modes = cm.get_cmap('tab10') # colormap for plotting mode network
    layer_color = {K:cmap_modes(i) for i,K in enumerate(K_range)}

    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)
    
    
    #%% Alignment within-K and mode detection
    
    tic = time.time()
    logging.info("---------- Aligning replicates within K and detecting modes ...")

    ILP_modes_filename = "ILPmodes.txt"
    
    if params.lc_flag:
        ## Leader Clustering
        if params.adaptive_thre_flag:

            logging.info(">>>Use leader clustering (adaptive threshold.")
            modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats = align_leader_clustering_adaptive(params.lc_cost_thre,Q_list,K_range,N,k2ids,ind2pop, pop_n_ind,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
        else:
            logging.info(">>>Use leader clustering (fixed threshold.")
            modes_allK_list,align_ILP_res,rep_modes,repQ_modes,meanQ_modes,average_stats = align_leader_clustering(params.lc_cost_thre,Q_list,K_range,N,k2ids,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
        
    else:
        logging.info(">>>Use ILP.")
        if params.cd_mod_thre!=-1:
            logging.info(">>>Use community detection modularity threshold: {})".format(params.cd_mod_thre))
        ## ILP over all pairs within K 
        # write to file then load from file
        ILP_withinK_filename = "ILPaligned.txt"
        align_ILP_withinK(ILP_withinK_filename,params.output_path,Q_list,K_range,k2ids)
        
        # load
        align_ILP_res, cost_ILP_res = load_ILP_withinK(Q_list,ILP_withinK_filename,output_path,K_range,k2ids,idx2idxinK,get_cost=True)
        
        ## Mode detection 
        modes_allK_list, msg = detect_modes(cost_ILP_res,K_range,params.default_cd,params.cd_mod_thre, save_path=save_path)
        if msg!="":
            logging.info(">>>"+msg)

        # extract representative/consensus replicates
        rep_modes,repQ_modes,meanQ_modes,average_stats = extract_modes(Q_list,modes_allK_list,align_ILP_res,cost_ILP_res,K_range,N,k2ids,save_modes=True, save_path=save_path,ILP_modes_filename=ILP_modes_filename)
    
    ## Within-K performance
    # stats of alignments: cost and hprime (similarity)
    weighted_stats = report_stats(modes_allK_list,average_stats,K_range,k2ids,save_path=save_path,stats=['cost','Hprime'])
    
    
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Alignment across-K
        
    # use representative membership as the consensus
    ILP_acrossK_filename = "ILP_acrossK_repQ.txt"
    logging.info("---------- Aligning modes (representative) across K ...")
    if not os.path.exists(os.path.join(save_path,ILP_acrossK_filename)):
        # write
        tic = time.time()
        repQ_acrossK_Q2P, repQ_acrossK_cost, repQ_best_ILP_acrossK = align_ILP_modes_acrossK(repQ_modes,K_range,N,save_path,ILP_acrossK_filename,enum_combK=params.enum_combK,ind2pop=ind2pop)
        toc = time.time()
        logging.info("Time: %.3fs", toc-tic)
    else:
        # load alignments 
        repQ_acrossK_Q2P, repQ_acrossK_cost, repQ_best_ILP_acrossK = load_ILP_acrossK(save_path,ILP_acrossK_filename)
        logging.info(">>>Across-K alignment (repQ) already exist.")
    

    # use average membership as the consensus
    ILP_acrossK_filename = "ILP_acrossK_meanQ.txt"
    logging.info("---------- Aligning modes (average) across K ...")
    if not os.path.exists(os.path.join(save_path+ILP_acrossK_filename)):
        # write
        tic = time.time()
        meanQ_acrossK_Q2P, meanQ_acrossK_cost, meanQ_best_ILP_acrossK = align_ILP_modes_acrossK(meanQ_modes,K_range,N,save_path,ILP_acrossK_filename,enum_combK=params.enum_combK,ind2pop=ind2pop)
        toc = time.time()
        logging.info("Time: %.3fs", toc-tic)
    else:
        meanQ_acrossK_Q2P, meanQ_acrossK_cost, meanQ_best_ILP_acrossK = load_ILP_acrossK(save_path,ILP_acrossK_filename)
        logging.info(">>>Across-K alignment (meanQ) already exist.")
    
    #%% Visualization of alignment within-K
    tic = time.time()
    logging.info("---------- Plotting alignment within K ...")
    ##  Alignment of each mode
    if params.plot_flag_all_within_mode:
        for K in K_range: 
            modes = modes_allK_list[K]
            for m in modes:
                plot_name = "K{}_mode{}_meanQ.pdf".format(K,m)
                plot_aligned(K,m,Q_list,modes,align_ILP_res,rep_modes,meanQ_modes[K][m],max_K,k2ids,idx2idxinK,save_path=save_path,plot_name=plot_name,cmap=cmap)
    
    ## Alignment between modes within-K
    if params.plot_flag_mode_within_K:
        plot_file_name_suffix = "modes"
        for K in K_range:
            plot_withinK_modes(K,max_K,meanQ_modes,meanQ_acrossK_Q2P,save_path,plot_file_name_suffix,cmap=cmap)
            
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Visualization of alignment across-K
    tic = time.time()
    logging.info("---------- Plotting alignment across K ...")
    ## Average membership multipartite graph
    if params.plot_flag_mode_across_K_multipartite:
        plot_file_name = "acrossK_meanQ_cost.pdf"
        title = "Mode (average membership) across K"
        G = plot_acrossK_multipartite(K_range,modes_allK_list,meanQ_modes,meanQ_acrossK_cost,layer_color,title,save_path,plot_file_name)
        
        # plot_file_name = "acrossK_repQ_cost.pdf"
        # title = "Mode (representative membership) across K"
        # G = plot_acrossK_multipartite(K_range,modes_allK_list,repQ_modes,repQ_acrossK_cost,layer_color,title,save_path,plot_file_name)
    
    ## Alignment of all modes
    if params.plot_flag_all_modes:
        plot_name = "all_modes_meanQ.pdf" 
        show_all_modes(params.plot_flag_all_modes,K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_best_ILP_acrossK,save_path,plot_name,cmap=cmap)    
        # plot_name = "all_modes_repQ.pdf" 
        # show_all_modes(params.plot_flag_all_modes,K_range,repQ_modes,repQ_acrossK_Q2P,repQ_best_ILP_acrossK,save_path,plot_name,cmap=cmap)    

    ## Optimal alignment across-K in chains
    if params.plot_flag_mode_across_K_chains:
        plot_file_name_suffix = "meanQ_alignment_chain_merged"
        best_alignment_chains = plot_acrossK_chains(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_acrossK_cost,save_path,plot_file_name_suffix,cmap=cmap,merge_cluster=True)
        
        plot_file_name_suffix = "meanQ_alignment_chain"
        best_alignment_chains = plot_acrossK_chains(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_acrossK_cost,save_path,plot_file_name_suffix,cmap=cmap,merge_cluster=False)
    
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    
    # zip all files
    logging.info("---------- Zipping files ...")
    shutil.make_archive(params.output_path, "zip", params.output_path)

    tot_toc = time.time()
    logging.info("======== Total Time: %.3fs ========", tot_toc-tot_tic)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input_path', type=str, required=True, help='path to the input files')
    required.add_argument('-o', '--output_path', type=str, required=True, help='path to the output files')
    required.add_argument('-p', '--params_path', type=str, required=True, help='path to the parameter file (.xml)')
    required.add_argument('-f', '--input_format', type=str, required=True, help='input data format')
    
    optional_arguments = [['default_cd','str','Y/N: whether to use default community detection method (Louvain) to detect modes'],
                          ['cd_mod_thre','float','the modularity threshold for community detection'],
                           ['custom_cmap','str','Y/N: whether to use customized colormap'],
                           ['cmap','str','user-specified colormap as a list of colors (in hex code) in a space-delimited string']]
    
    for opt_arg in optional_arguments: 
        optional.add_argument('--{}'.format(opt_arg[0]), type=getattr(builtins, opt_arg[1]), required=False, help=opt_arg[2])

    args = parser.parse_args()
    main(args)