# -*- coding: utf-8 -*-
"""
Clumppling diffInd function: not yet ready for beta version.

@author: Xiran Liu
"""

import numpy as np
import os
import sys
import matplotlib.cm as cm
import time
import logging
import shutil
import builtins
import matplotlib.pyplot as plt
import argparse

from clumppling.funcs import load_parameters, get_random_cmap, load_Q_and_indinfo, align_popwise_membership


def main(args):
    
    input_base_path = args.input_path
    output_path = args.output_path

    # sanity check for arguments
    if not os.path.exists(input_base_path):
        sys.exit("ERROR: Input file path doesn't exist.")

    # parameters = load_parameters(params_path)
    
    # create output directory
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    output_f = os.path.join(output_path,'output.log')
    handlers = [logging.FileHandler(output_f, 'w'), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format='', handlers = handlers)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True
    
    logging.info("========== [Parameters] ========== ")
    logging.info("Input path: {}".format(input_base_path))
    logging.info("Output path: {}".format(output_path))
    logging.info("Parameter file path: {}".format(params_path))
    
    tot_tic = time.time()

    logging.info("========== Running Clumppling (diffInd) ==========")
    
    # determine input directories
    input_names = [f for f in os.listdir(input_base_path) if os.path.isdir(os.path.join(input_base_path, f))]
    logging.info("Input Q files with different individuals: {}".format(", ".join(input_names)))
    logging.info("Consensus: {}".format(parameters['cons_suffix']))

    # load files
    tic = time.time()
    logging.info("---------- Loading files ...")
    input_names, N_all, R_all, Q_all, K_all, ind2pop_all, popNind_all = load_Q_and_indinfo(input_base_path,parameters['cons_suffix'])
    K_max = max(max(k) for k in K_all)
    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)

    # set colormap
    if parameters['custom_cmap']!='':
        custom_cmap = parameters['custom_cmap'].split(",")
        while len(custom_cmap) < K_max:
            logging.info(">>>The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            custom_cmap.extend(custom_cmap)
        cmap = cm.colors.ListedColormap(custom_cmap)
    else:
        cmap = get_random_cmap(cm.get_cmap('Spectral'),K_max,seed=999)

    
    # alignment
    tic = time.time()
    logging.info("---------- Aligning and plotting ...")


    # align and plot
    align_popwise_membership(input_names, Q_all, K_all, ind2pop_all, output_path, cmap)

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

    required.add_argument('-i', '--input_path', type=str, required=True, help='path to the input files')
    required.add_argument('-o', '--output_path', type=str, required=True, help='path to the output files')
    
    
    args = parser.parse_args()
    main(args)