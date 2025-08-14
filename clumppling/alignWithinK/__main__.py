import os
import argparse
from typing import List
from clumppling.log_config import setup_logger
from clumppling.utils import disp_params, disp_msg, parse_strings
from clumppling.core import align_within_k
from . import load_withinK_qfiles, write_aligned_within_k


def parse_args():
    parser = argparse.ArgumentParser(description="clumppling.alignWithinK")

    parser.add_argument("--qfiles", nargs='*', help="List of Q files to align, passed as command-line arguments")
    parser.add_argument(
        "--qfilelist", type=str,
        help="A plain text file containing Q file names (one per line).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file name")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(os.path.dirname(args.output), "alignWithinK.log"))
    disp_params(args, title="Alignment within K")

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    disp_msg(f"Loading Q files for alignment within K")
    Q_files = parse_strings(args.qfilelist, args.qfiles) #parse_qfiles(args)
    Q_names = [os.path.splitext(os.path.basename(p))[0] for p in Q_files]
    Q_list = load_withinK_qfiles(Q_files)  # Example K value, replace with actual logic to determine K
    disp_msg(f"Aligning {len(Q_list)} Q files within K")
    alined_results = align_within_k(Q_list)
    disp_msg(f"Writing aligned results to {args.output}")
    write_aligned_within_k(alined_results, Q_names, args.output)
