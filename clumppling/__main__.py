# -*- coding: utf-8 -*-
"""
Clumppling main function

@author: Xiran Liu
"""

import numpy as np
import os
import sys
import json
# import matplotlib.cm as cm
import time
import logging
import shutil
import builtins   
import argparse
from pkg_resources import resource_stream

from clumppling.funcs import *
# import warnings
# warnings.filterwarnings("ignore")



def main(args):

    input_path = args.input_path
    output_path = args.output_path
    # params_path = args.params
    input_format = args.input_format
    
    # sanity check for arguments
    if os.path.exists(output_path):
        sys.exit("ERROR: Output directory {} already exists. Please remove the directory before running the program.".format(output_path))
    if os.path.exists(input_path+".zip"):
        shutil.unpack_archive(input_path+".zip",input_path)
    if not os.path.exists(input_path):
        sys.exit("ERROR: Input file doesn't exist.")
    if not input_format in ["structure","fastStructure","admixture","generalQ"]:
        sys.exit("ERROR: Input data format is not supported. \nPlease specify input_format as one of the following: structure, admixture, fastStructure, and generalQ.")

    if args.vis is not None:
        visualization = bool(args.vis)
    else:
        visualization = True


    if args.params is None:
        params_path = 'default parameters'
        params_default = json.load(resource_stream('clumppling', 'files/default_params.json'))
        parameters = load_default_parameters(params_default)
    else:
        params_path = args.params_path
        if not os.path.exists(params_path):
            sys.exit("ERROR: Parameter file doesn't exist.")
        parameters = load_parameters(params_path)
    
    if args.cd_param is not None:
        parameters['cd_parameter'] = args.cd_param if args.cd_param>0 else 1.0
    if args.use_rep is not None:
        parameters['use_rep'] = bool(args.use_rep) 
    if args.merge_cls is not None:
        parameters['merge_cls'] = bool(args.merge_cls) 
    plot_params = ['plot_modes','plot_modes_withinK','plot_major_modes','plot_all_modes']
    if visualization==False:
        for param in plot_params:
            parameters[param] = False
    disp = display_parameters(input_path,input_format,output_path,params_path,parameters)
    
    #%% Set-up 
    tot_tic = time.time()
    
    # if os.path.exists(output_path):
    #     shutil.rmtree(output_path)
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
    logging.info("======= Running Clumppling =======")


    tic = time.time()
    logging.info(">>> Processing input data files and checking arguments")
    
    
    #%% Loading membership data    
    # load input data
    Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs = load_inputs(data_path=input_path, 
                                                                   output_path=output_path, 
                                                                   input_format=input_format)
    
    distruct_cmap = [
        '#FF994D','#0099E6','#E6FF00','#FF99E6','#339933','#800080','#FF004D','#00FF00','#0000FF','#FF00FF',
        '#FFE699','#B24D00','#00FFFF','#808000','#FF9999','#008080','#99BF26','#7326E6','#26BF99','#808080',
        '#0D660D','#BFBFBF','#FF0000','#99E6FF','#FF9966','#404040','#FFE6E6','#993333','#FF6600','#33004D',
        '#FFFFFF','#FF4D00','#FF9900','#FF4D99','#FF99E6','#FFFF99','#FFFFE6','#FEB1B1','#E6FF4D','#E6E6FF',
        '#E699FF','#99FF4D','#99FF99','#99FFE6','#99FFFF','#99E6FF','#9999FF','#994DFF','#4D99FF','#00FF99',
        '#00FFE6','#00E6FF','#CC1A1A','#B24D00','#B2331A','#B21A33','#996600','#994D1A','#4D664D','#4D664D']
    
    # set colormap
    if parameters['custom_cmap']!='':
        custom_cmap = [s.strip() for s in parameters['custom_cmap'].split(",")]
        while len(custom_cmap) < K_max:
            logging.info("The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            custom_cmap.extend(custom_cmap)
        cmap = cm.colors.ListedColormap(custom_cmap)
    else:
        cmap = cm.colors.ListedColormap(distruct_cmap[:K_max])

    if visualization:
        # visualize results
        fig_path = os.path.join(output_path,"visualization")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        plot_colorbar(cmap,K_max,fig_path)


    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)
    
    
    #%% Alignment within-K and mode detection
    
    tic = time.time()
    logging.info(">>> Aligning replicates within K and detecting modes")

    # align within-K
    alignment_withinK, cost_withinK = align_withinK(output_path,Q_list,Q_files,K_range,K2IDs)

    # detect mode
    modes_allK, cost_matrices, msg = detect_modes(cost_withinK,Q_files,K_range,K2IDs,default_cd=parameters['cd_default'],cd_mod_thre=parameters['cd_modularity_threshold'],cd_param=parameters['cd_parameter'])
    
    # if msg!="":
    #     logging.info("Mode detection finishes.\n"+msg)

    # extract representative/consensus replicates
    mode_labels, rep_modes, repQ_modes, avgQ_modes, alignment_to_modes, stats = extract_modes(Q_list,Q_files,modes_allK,alignment_withinK,cost_matrices,K_range,K2IDs, output_path)
    
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Alignment across-K
        
    # use representative membership as the consensus
    logging.info(">>>  Aligning modes across K")
    tic = time.time()
    # align across-K
    acrossK_path = os.path.join(output_path,"alignment_acrossK")
    alignment_acrossK_avg, cost_acrossK_avg, best_acrossK_avg = align_ILP_modes_acrossK(avgQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="avg",merge=parameters['merge_cls'])
    alignment_acrossK_rep, cost_acrossK_rep, best_acrossK_rep = align_ILP_modes_acrossK(repQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="rep",merge=parameters['merge_cls'])
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Visualization of alignment within-K
    tic = time.time()
    logging.info(">>> Plotting alignment results")
    if visualization and parameters['plot_modes_withinK']:
        for K in K_range:
            if parameters['use_rep']:
                plot_withinK_modes(K,K_max,repQ_modes,alignment_acrossK_rep,fig_path,cmap,fig_suffix="rep")
            else:
                plot_withinK_modes(K,K_max,avgQ_modes,alignment_acrossK_avg,fig_path,cmap,fig_suffix="avg")
    

    #%% Visualization of alignment across-K
    # # plot all replicates
    # plot_replicates(Q_list,K_range,output_path,cmap)
    
    if visualization and parameters['plot_modes']:
        if parameters['use_rep']:
            plot_structure_on_multipartite(K_range,mode_labels,stats,repQ_modes,alignment_acrossK_rep,cost_acrossK_rep,best_acrossK_rep,"rep",output_path,True,cmap)
        else:
            # plot_structure_on_multipartite_manuscript(K_range,mode_labels,stats,avgQ_modes,alignment_acrossK_avg,cost_acrossK_avg,best_acrossK_avg,"avg",output_path,True,cmap)
            plot_structure_on_multipartite(K_range,mode_labels,stats,avgQ_modes,alignment_acrossK_avg,cost_acrossK_avg,best_acrossK_avg,"avg",output_path,True,cmap)
    else:
        if parameters['use_rep']:
            plot_structure_on_multipartite(K_range,mode_labels,stats,repQ_modes,alignment_acrossK_rep,cost_acrossK_rep,best_acrossK_rep,"rep",output_path,False,cmap)
        else:
            plot_structure_on_multipartite(K_range,mode_labels,stats,avgQ_modes,alignment_acrossK_avg,cost_acrossK_avg,best_acrossK_avg,"avg",output_path,False,cmap)

    if visualization and parameters['plot_major_modes']:
        if parameters['use_rep']:
            plot_major_modes(K_range,output_path,"rep",cmap)
        else:
            plot_major_modes(K_range,output_path,"avg",cmap)
    if visualization and parameters['plot_all_modes']:
        if parameters['use_rep']:
            plot_all_modes(K_range,mode_labels,output_path,"rep",cmap)
        else:
            plot_all_modes(K_range,mode_labels,output_path,"avg",cmap)

    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    
    # zip all files
    logging.info(">>>  Zipping files")
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
    required.add_argument('-f', '--input_format', type=str, required=True, help='input data format')
    optional.add_argument('-p', '--params', type=str, required=False, help='path to the parameter file (.json)')
    
    
    optional.add_argument('-v', '--vis', type=int, required=False, help='whether to generate visualization: 0 for no, 1 for yes (default)')
    
    optional_arguments = [['cd_param','float','the parameter for community detection method'], \
        ['use_rep','int','whether to use representative replicate as mode consensus: 0 for no (default), 1 for yes'], \
        ['merge_cls','int','whether to merge all pairs of clusters to align K+1 and K: 0 for no (default), 1 for yes']]
    for opt_arg in optional_arguments: 
        optional.add_argument('--{}'.format(opt_arg[0]), type=getattr(builtins, opt_arg[1]), required=False, help=opt_arg[2])

    args = parser.parse_args()

    main(args)