import os
import argparse
from clumppling.log_config import setup_logger
from clumppling.utils import disp_msg, str2bool, disp_params
from . import process_files

def parse_args():
    parser = argparse.ArgumentParser(description="clumppling.parseInput")

    # required
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file path")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file directory")
    parser.add_argument("-f", "--format", type=str, required=True, choices=["generalQ", "admixture", "structure", "fastStructure"],
                        help="File format")
    # optional
    parser.add_argument("--extension", type=str, default="", help="Extension of input files")
    parser.add_argument("--skip_rows", type=int, default=0, help="Skip top rows in input files")
    parser.add_argument("--remove_missing", type=str2bool, default=True, help="Remove individuals with missing data: True/False")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(args.output, "parseInput.log"))
    disp_params(args, title="Input Parsing")
    disp_msg(f"Processing input files")
    process_files(input_dir=args.input, output_dir=os.path.join(args.output,"input"), 
                  fmt=args.format,extension=args.extension, 
                  skip_missing=args.remove_missing,
                  delimiter=" ", skip_rows=0)
