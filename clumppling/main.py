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
import argparse

from clumppling.funcs import *



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
        sys.exit("ERROR: Input data format is not supported. \nPlease specify input_format as one of the following: structure, admixture, fastStructure, and generalQ.")
    
    parameters = load_parameters(params_path)
    if args.cd_mod_thre:
        parameters['cd_modularity_threshold'] = args.cd_mod_thre if args.cd_mod_thre!=-1 else -1.0
    disp = display_parameters(input_path,input_format,output_path,params_path,parameters)
    
    
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
    
    output_f = os.path.join(output_path,'output.log')
    handlers = [logging.FileHandler(output_f, 'w'), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format='', handlers = handlers)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True

    logging.info(disp)
    logging.info("========== Running Clumppling ==========")


    tic = time.time()
    logging.info("---------- Processing input data files and checking arguments ...")
    
    
    #%% Loading membership data    
    # load input data
    Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs = load_inputs(data_path=input_path, 
                                                                   output_path=output_path, 
                                                                   input_format=input_format)
    
    # set colormap
    if parameters['custom_cmap']!='':
        custom_cmap = parameters['custom_cmap'].split(",")
        while len(custom_cmap) < K_max:
            logging.info(">>>The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            custom_cmap.extend(custom_cmap)
        cmap = cm.colors.ListedColormap(custom_cmap)
    else:
        cmap = get_random_cmap(cm.get_cmap('Spectral'),K_max,seed=999)

    # visualize results
    fig_path = os.path.join(output_path,"visualization")
    if not os.path.exists(fig_path):
        os.mkdir(fig_path)

    plot_colorbar(cmap,K_max,fig_path)


    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)
    
    
    #%% Alignment within-K and mode detection
    
    tic = time.time()
    logging.info("---------- Aligning replicates within K and detecting modes ...")

    # align within-K
    alignment_withinK, cost_withinK = align_withinK(output_path,Q_list,Q_files,K_range,K2IDs)

    # detect mode
    modes_allK, cost_matrices, msg = detect_modes(cost_withinK,Q_files,K_range,K2IDs,cd_default,cd_mod_thre=parameters['cd_modularity_threshold'],cd_param=parameters['cd_parameter'])
    
    if msg!="":
        logging.info(">>>"+msg)

    # extract representative/consensus replicates
    mode_labels, rep_modes, repQ_modes, meanQ_modes, alignment_to_modes, stats  = extract_modes(Q_list,Q_files,modes_allK,alignment_withinK,cost_matrices,K_range,K2IDs, output_path)
    
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Alignment across-K
        
    # use representative membership as the consensus
    logging.info("---------- Aligning modes across K ...")
    tic = time.time()
    # align across-K
    acrossK_path = os.path.join(output_path,"alignment_acrossK")
    alignment_acrossK_mean, cost_acrossK_mean, best_acrossK_mean = align_ILP_modes_acrossK(meanQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="mean",merge=False)
    alignment_acrossK_rep, cost_acrossK_rep, best_acrossK_rep = align_ILP_modes_acrossK(repQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="rep",merge=False)
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Visualization of alignment within-K
    tic = time.time()
    logging.info("---------- Plotting alignment within K ...")
    if parameters['plot_modes_withinK']:
        for K in K_range:
            plot_withinK_modes(K,K_max,meanQ_modes,alignment_acrossK_mean,fig_path,cmap,fig_suffix="mean")
            plot_withinK_modes(K,K_max,repQ_modes,alignment_acrossK_rep,fig_path,cmap,fig_suffix="rep")
      
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Visualization of alignment across-K
    tic = time.time()
    logging.info("---------- Plotting alignment across K ...")
    if parameters['plot_modes']:
        align_all_modes(K_range,mode_labels,meanQ_modes,alignment_acrossK_mean,best_acrossK_mean,output_path,"mean",True,cmap)
        align_all_modes(K_range,mode_labels,repQ_modes,alignment_acrossK_rep,best_acrossK_rep,output_path,"rep",True,cmap)
    else:
        align_all_modes(K_range,mode_labels,meanQ_modes,alignment_acrossK_mean,best_acrossK_mean,output_path,"mean",False,cmap)
        align_all_modes(K_range,mode_labels,repQ_modes,alignment_acrossK_rep,best_acrossK_rep,output_path,"rep",False,cmap)
    
    if parameters['plot_acrossK']:
        plot_acrossK_multipartite(K_range,mode_labels,stats,cost_acrossK_mean,fig_path,"mean")
        plot_acrossK_multipartite(K_range,mode_labels,stats,cost_acrossK_rep,fig_path,"rep")

    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    
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

    required.add_argument('-i', '--input_path', type=str, required=True, help='path to the input files')
    required.add_argument('-o', '--output_path', type=str, required=True, help='path to the output files')
    required.add_argument('-p', '--params_path', type=str, required=True, help='path to the parameter file (.json)')
    required.add_argument('-f', '--input_format', type=str, required=True, help='input data format')
    
    optional_arguments = [['cd_mod_thre','float','the modularity threshold for community detection (default: -1, not using the threshold)']]

    for opt_arg in optional_arguments: 
        optional.add_argument('--{}'.format(opt_arg[0]), type=getattr(builtins, opt_arg[1]), required=False, help=opt_arg[2])

    args = parser.parse_args()
    main(args)