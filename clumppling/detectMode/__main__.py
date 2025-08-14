import os
import pandas as pd
import argparse
import shutil
from typing import List
from . import construct_cost_mat, write_modes_to_file, community_labels_to_modes

from clumppling.log_config import setup_logger
from clumppling.utils import disp_params, disp_msg, str2bool, parse_strings
from clumppling.core import detect_communities, extract_modes_and_stats
from clumppling.alignWithinK import load_aligned_within_k, aligned_res_to_dicts, load_withinK_qfiles

import logging
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="clumppling.alignWithinK")

    parser.add_argument("--align_res", type=str, required=True,
                        help="Path to the alignment results file")
    parser.add_argument("--qfilelist", type=str, required=True,
                        help="A plain text file containing Q file names (one per line).")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Output file directory")
    parser.add_argument("--qnamelist", type=str, default=[],
                        help="A plain text file containing replicate names (one per line) (default: file base from qfilelist)")
    parser.add_argument("--cd_method", type=str, default="louvain",
                        choices=["louvain", "leiden", "infomap", "markov_clustering", "label_propagation", "walktrap", "custom"],
                        help="Community detection method to use (default: louvain)")
    parser.add_argument("--cd_res", type=float, default=1.0,
                        help="Resolution parameter for the default Louvain community detection (default: 1.0)")
    parser.add_argument("--test_comm", type=str2bool, default=True,
                        help="Whether to test community structure (default: True)")
    parser.add_argument("--comm_min", type=float, default=1e-4,
                        help="Minimum threshold for cost matrix (default: 1e-4)")
    parser.add_argument("--comm_max", type=float, default=1e-2,
                        help="Maximum threshold for cost matrix (default: 1e-2)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(args.output, "detectMode.log"))
    disp_params(args, title="Mode Detection")

    if os.path.exists(args.output) and os.listdir(args.output):
        shutil.rmtree(args.output)
        logger.info(f"Output directory '{args.output}' already exists and is not empty. Removed existing directory.")
    os.makedirs(args.output, exist_ok=True)
    logger.info(f"Created output directory '{args.output}'.")

    disp_msg(f"Loading replicate names")
    Q_files = parse_strings(args.qfilelist) 
    if args.qnamelist is None:
        Q_names = [os.path.splitext(os.path.basename(p))[0] for p in Q_files]
    else:
        Q_names = parse_strings(args.qnamelist)
    disp_msg(f"Loading alignment within-K results")
    aligned_summary = load_aligned_within_k(input_file = args.align_res, Q_names = Q_names)
    alignment_withinK, cost_withinK = aligned_res_to_dicts(aligned_summary, Q_names)
    disp_msg(f"Constructing cost matrix")
    cost_mat = construct_cost_mat(cost_withinK, Q_names)   
    disp_msg(f"Detecting communities")
    communities = detect_communities(cost_mat, test_comm = args.test_comm,
                                     method = args.cd_method, res = args.cd_res,
                                     min_threshold = args.comm_min, max_threshold = args.comm_max)   
    modes = community_labels_to_modes(communities)
    disp_msg(f"Extracting modes (Q matrices) and statistics")
    Q_list = load_withinK_qfiles(Q_files)
    res = extract_modes_and_stats(modes, cost_mat, Q_names, 
                                  Q_list, alignment_withinK)
    disp_msg(f"Saving modes and statistics")
    avg_stat = write_modes_to_file(res, output_dir=args.output, compute_avg=True)
    print(avg_stat)

    
    