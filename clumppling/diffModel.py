# -*- coding: utf-8 -*-
"""
Clumppling diffModel function: not yet ready for beta version.

@author: Xiran Liu
"""

import os
import sys
import time
import shutil
import logging
import datetime
import builtins
import argparse
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pkg_resources import resource_stream


from clumppling.funcs import load_Q_and_indinfo, align_multiple_model


def main(args):
    
    input_base_path = args.input_path
    output_path = args.output_path

    # sanity check for arguments
    if not os.path.exists(input_base_path):
        sys.exit("ERROR: Input file path {} doesn't exist.".format(input_base_path))

    parameters = vars(args)    
    if not parameters['consensus'] in ["avg","rep"]:
        sys.exit("ERROR: mode consensus is not recognized. \nPlease specify 'consensus' as either 'avg' or 'rep'.")  

    # create output directory
    if os.path.exists(output_path+".zip"):
        os.remove(output_path+".zip")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    #%% Set-up 
    # log outputs
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=os.path.join(output_path,'output.log'), level=logging.INFO, format='')
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger('matplotlib.pyplot').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True

    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info("{} Program starts.".format(current_time))

    logging.info("==================================")
    logging.info("Input path: {}".format(input_base_path))
    logging.info("Output path: {}".format(output_path))
    logging.info("= Running Clumppling (diffModel) =")

    tot_tic = time.time()

    # determine input directories
    input_names = [f for f in os.listdir(input_base_path) if os.path.isdir(os.path.join(input_base_path, f))]
    logging.info("Input Q files with different models: {}".format(", ".join(input_names)))
    logging.info("Mode consensus: {}".format(parameters['consensus']))

    # load files
    tic = time.time()
    logging.info("---------- Loading files ...")
    input_names, N_all, R_all, Q_all, K_all = load_Q_and_indinfo(input_base_path,parameters['consensus'],indinfo=False)
    K_max = max(max(k) for k in K_all)
    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)

    # set colormap
    distruct_cmap = np.loadtxt(resource_stream('clumppling', 'files/default_colors.txt'),dtype=str,comments=None)
    if parameters['custom_cmap']!='':
        custom_cmap = [s.strip() for s in parameters['custom_cmap'].split(",")]
        while len(custom_cmap) < K_max:
            logging.info("The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            custom_cmap.extend(custom_cmap)
        cmap = cm.colors.ListedColormap(custom_cmap)
    else:
        cmap = cm.colors.ListedColormap(distruct_cmap[:K_max])

    
    # alignment
    tic = time.time()
    logging.info("---------- Aligning and plotting ...")

    # align and plot
    align_multiple_model(input_names, Q_all, K_all, output_path, cmap)
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

    required.add_argument('-i', '--input_path', type=str, required=True, help='path to the input files')
    required.add_argument('-o', '--output_path', type=str, required=True, help='path to the output files')
    optional.add_argument('--consensus', default='avg', type=str, required=False, help='what to use as mode consensus, either avg (default, for average memeberships) or rep (for representative replicate)')
    optional.add_argument('--custom_cmap', default='', type=str, required=False, help='customized colormap as a comma-separated string of hex codes for colors: if empty (default), using the default colormap, otherwise use the user-specified colormap')

    args = parser.parse_args()
    main(args)