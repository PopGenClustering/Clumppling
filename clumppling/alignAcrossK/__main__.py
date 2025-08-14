import os
import shutil
import argparse

from clumppling.log_config import setup_logger
from clumppling.utils import disp_params, disp_msg, parse_strings
from clumppling.core import align_across_k, write_alignment_across_k
from . import load_any_qfiles, extract_K_range_from_Qs, separate_Qs_by_K

import logging
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="clumppling.alignAcrossK")

    parser.add_argument(
        "--qfilelist", type=str,
        help="A plain text file containing Q file names (one per line).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Directory to save output files")
    parser.add_argument("--qnamelist", type=str, default=[],
                        help="A plain text file containing replicate names (one per line) (default: file base from qfilelist)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(args.output, "alignAcrossK.log"))
    disp_params(args, title="Alignment across K")

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
    disp_msg(f"Extracting K values from Q files")
    Q_list = load_any_qfiles(Q_files)
    K_range = extract_K_range_from_Qs(Q_list)
    if len(K_range) < 2:
        logger.error("At least two K values are required for alignment across K.")
    Q_modes_list, mode_names_list = separate_Qs_by_K(K_range, Q_list, Q_names)
    disp_msg(f"Aligning")
    alignment_acrossK, cost_acrossK, best_acrossK_out, major_acrossK_out = align_across_k(K_range, Q_modes_list, mode_names_list, merge=True)
    disp_msg(f"Saving alignment results to {args.output}")
    write_alignment_across_k(alignment_acrossK, cost_acrossK, os.path.join(args.output,"alignment_acrossK.txt"))
    best_acrossK_out.to_csv(os.path.join(args.output,"best_pairs_acrossK.txt"), index=False)
    major_acrossK_out.to_csv(os.path.join(args.output,"major_pairs_acrossK.txt"), index=False)
    

