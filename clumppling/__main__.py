import argparse
import os
import shutil
import time
import matplotlib.pyplot as plt

from .log_config import setup_logger
from .core import align_within_k_all_K, detect_modes_all_K, extract_modes_all_K, align_across_k, write_alignment_across_k, reorderQ_within_k, reorderQ_across_k
from .utils import disp_params,str2bool,parse_strings,get_modes_all_K, write_reordered_across_k, get_mode_sizes, unnest_list
from .plot import load_default_cmap, parse_custom_cmap, plot_colorbar, plot_memberships_list, plot_graph, plot_alignment

from .parseInput import process_files, extract_meta_input, group_labels

import logging
logger = logging.getLogger(__name__)

def main(args: argparse.Namespace):
    tot_tic = time.time()

    # load and process input files
    processed_input_dir = os.path.join(args.output, "input")
    if os.path.exists(processed_input_dir) and os.listdir(processed_input_dir):
        shutil.rmtree(processed_input_dir)
        logger.info(f"Directory '{processed_input_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(processed_input_dir, exist_ok=True)
    logger.info(f"Created input directory '{processed_input_dir}'.")

    # load and reorder individuals labels if applicable
    reorder_indices = None
    if args.include_label:
        if args.ind_labels!="":
            ind_labels = parse_strings(args.ind_labels, remove_dup=False)
            if args.regroup_ind and len(ind_labels)>0:
                logger.info(f"{len(ind_labels)} individual labels loaded.")
                ind_labels, reorder_indices = group_labels(ind_labels)
                if not reorder_indices == list(range(len(ind_labels))):
                    logger.warning(f"Individual labels reordered by grouping identical labels together.")
                    # Save grouped labels
                    with open(os.path.join(processed_input_dir,"ind_labels_grouped.txt"), "w") as f:
                        for label in ind_labels:
                            f.write(label + "\n")
                    # Save reorder indices (0-based)
                    with open(os.path.join(processed_input_dir,"ind_labels_indices.txt"), "w") as f:
                        for i in reorder_indices:
                            f.write(str(i) + "\n")
            else:
                logger.info(f"{len(ind_labels)} individual labels loaded (not regrouped).")

    # process files
    labels = process_files(input_dir=args.input, output_dir=processed_input_dir, 
                           fmt=args.format,extension=args.extension, 
                           skip_missing=args.remove_missing,
                           delimiter=" ", skip_rows=0, label_cols=[0, 1, 3], 
                           mat_start_col=5, reorder_indices=reorder_indices)

    if args.include_label:
        if args.ind_labels!="":
            logger.info(f"Input labels loaded from {args.ind_labels}.")
        elif labels is not None:
            ind_labels = list(labels[:,-1])
            logger.info("Input labels loaded from processed files.")
        else:
            ind_labels = []
            logger.warning("No input labels found.")
    else:
        ind_labels = []
        logger.warning("No input labels.")
    

    Q_names, K_range, K2IDs = extract_meta_input(processed_input_dir)
    K_max = max(K_range)
    logging.info(f"Unique K values found (max: {K_max}): {K_range}")

    # align within K
    logger.info(f"Aligning replicates within K".center(50, '-'))
    alignment_withinK_list, cost_withinK_list = align_within_k_all_K(Q_names, K_range, K2IDs,
                                                                     qfile_dir = processed_input_dir, output_dir = args.output)
    # detect modes within K
    logger.info(f"Detecting modes within K".center(50, '-'))
    modes_dir = os.path.join(args.output, "modes")
    modes_all_K_list, cost_matrices_list = detect_modes_all_K(K_range, cost_withinK_list, Q_names, K2IDs, 
                                                    test_comm = args.test_comm, method = args.cd_method, res = args.cd_res, 
                                                    comm_min = args.comm_min, comm_max = args.comm_max)
    # extract modes and statistics
    logger.info(f"Extracting modes and summarizing statistics".center(50, '-'))
    cd_res, avg_stat = extract_modes_all_K(K_range, K2IDs, Q_names, cost_matrices_list, modes_all_K_list, alignment_withinK_list,
                                           processed_input_dir, output_dir=modes_dir)
    mode_sizes = get_mode_sizes(cd_res) 
    # get modes for alignment across K
    mode_names_list, Q_rep_modes_list, Q_avg_modes_list = get_modes_all_K(K_range, cd_res)
    
    # align across K
    logger.info(f"Aligning modes across K".center(50, '-'))
    if len(K_range) < 2:
        logger.error("At least two K values are required for alignment across K: K={}".format(K_range))
    acrossK_dir = os.path.join(args.output,"alignment_acrossK")
    if os.path.exists(acrossK_dir) and os.listdir(acrossK_dir):
        shutil.rmtree(acrossK_dir)
        logger.info(f"Alignment across-K output directory '{acrossK_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(acrossK_dir, exist_ok=True)
    suffix = "rep" if args.use_rep else "avg"
    if args.use_rep:
        logger.info(f"(using representative replicate)")
        alignment_acrossK, cost_acrossK, best_acrossK_out, major_acrossK_out = align_across_k(K_range, Q_rep_modes_list, mode_names_list, merge=args.merge)
        write_alignment_across_k(alignment_acrossK, cost_acrossK, os.path.join(acrossK_dir,"alignment_acrossK_rep.txt"))
        best_acrossK_out.to_csv(os.path.join(acrossK_dir,"best_pairs_acrossK_rep.txt"), index=False)
        major_acrossK_out.to_csv(os.path.join(acrossK_dir,"major_pairs_acrossK_rep.txt"), index=False)
    else:
        logger.info(f"(using average across replicates)")
        alignment_acrossK, cost_acrossK, best_acrossK_out, major_acrossK_out  = align_across_k(K_range, Q_avg_modes_list, mode_names_list, merge=args.merge)
        write_alignment_across_k(alignment_acrossK, cost_acrossK, os.path.join(acrossK_dir,"alignment_acrossK_rep.txt"))
        best_acrossK_out.to_csv(os.path.join(acrossK_dir,"best_pairs_acrossK_avg.txt"), index=False)
        major_acrossK_out.to_csv(os.path.join(acrossK_dir,"major_pairs_acrossK_avg.txt"), index=False)    
    
    # Reorder Q matrices across K
    logger.info(f"Reordering membership matrices according to alignment".center(50, '-'))
    modeQ_dir = os.path.join(args.output,"modes_aligned")
    Q_modes_list = Q_rep_modes_list if args.use_rep else Q_avg_modes_list
    anchor_pairs = best_acrossK_out["Best Pair"].tolist() if args.use_best_pair else major_acrossK_out["Major Pair"].tolist()
    aligned_Qs_allK, all_modes_alignment = reorderQ_across_k(K_range, Q_modes_list, mode_names_list, 
                                                            alignment_acrossK, anchor_pairs)
    write_reordered_across_k(aligned_Qs_allK, all_modes_alignment, output_dir=modeQ_dir, suffix=suffix)
    if len(K_range)==1:
        logger.warning(f"Only one K value found: K={K_range[0]}." )
    
    # plot colormap if visualization is enabled
    if args.vis:
        # load colormap (convert to a list of RGB tuples)
        if args.custom_cmap:
            try:
                colors = parse_strings(args.custom_cmap, remove_dup=False)
                cmap = parse_custom_cmap(colors, K=K_max)
            except ValueError as e:
                logger.error(e)
                logger.info("Falling back to default colormap.")
                cmap = load_default_cmap(K=K_max)
        else:
            cmap = load_default_cmap(K=K_max)  # Example usage of load_default_cmap
            
        # plot alignment pattern
        logger.info(f"Generating visualizations".center(50, '-'))
        fig_dir = os.path.join(args.output,"visualization")
        os.makedirs(fig_dir, exist_ok=True)
        logger.info(f"Plot colorbar")
        plot_colorbar(cmap,K_max,fig_dir)

        # plot alignment pattern
        logger.info(f"Plot alignment patterns ({suffix})")
        mode_K = [Q.shape[1] for Q in unnest_list(Q_rep_modes_list)]
        mode_names = unnest_list(mode_names_list)
        fig = plot_alignment(mode_K, mode_names, cmap, alignment_acrossK, all_modes_alignment, marker_size=200)
        fig.savefig(os.path.join(fig_dir,"alignment_pattern_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig)

        # determine width scaling for structure plots based on number of unique individual labels
        width_scale = 1.0
        if len(ind_labels)>1:
            width_scale = max(1,len(set(ind_labels))/8)

        if len(K_range)>1:

            if args.plot_type in ["graph", "all"]:
                logger.info(f"Plot all modes in a graph ({suffix})")
                mode_labels_list = [[f"{mode_name} ({mode_sizes[mode_name]})" for mode_name in mode_names] for mode_names in mode_names_list]
                Q_modes_reordered_list = [[aligned_Qs_allK[mode_name] for mode_name in mode_names] for mode_names in mode_names_list]
                if args.reorder_ind:
                    if args.reorder_by_max_k:
                        Q_ref = Q_modes_reordered_list[-1][0]
                    else:
                        Q_ref = Q_modes_reordered_list[0][0]
                else:
                    Q_ref = None
                
                if args.include_cost:
                    logger.info(f"Including cost values in the graph")
                    fig = plot_graph(K_range, Q_modes_reordered_list, cmap, 
                                    names_list=mode_names_list, labels_list=mode_labels_list,  
                                    cost_acrossK=cost_acrossK, ind_labels=ind_labels, 
                                    fontsize=14, line_cmap=plt.get_cmap("Greys"),
                                    order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label,
                                    width_scale=width_scale)
                else:
                    logger.info(f"Not including cost values in the graph")
                    fig = plot_graph(K_range, Q_modes_reordered_list, cmap, 
                                    names_list=mode_names_list, labels_list=mode_labels_list,  
                                    cost_acrossK=None, ind_labels=ind_labels, 
                                    fontsize=14, line_cmap=plt.get_cmap("Greys"),
                                    order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label,
                                    width_scale=width_scale)
                fig.savefig(os.path.join(fig_dir,"all_modes_graph_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
                plt.close(fig)
            
            if args.plot_type in ["withinK", "all"]:
                logger.info(f"Plot modes within-K ({suffix})")
                for i_K, K in enumerate(K_range):
                    Q_modes_list = Q_rep_modes_list[i_K] if args.use_rep else Q_avg_modes_list[i_K]
                    mode_names = mode_names_list[i_K]
                    Q_modes_reordered = [aligned_Qs_allK[mode_name] for mode_name in mode_names]
                    # # if only align within each mode, but not across K, then use the following line:
                    # Q_modes_reordered = reorderQ_within_k(Q_modes_list, mode_names, alignment_acrossK)
                    mode_labels = [f"{mode_name} ({mode_sizes[mode_name]})" for mode_name in mode_names]
                    if args.reorder_ind:
                        if args.reorder_by_max_k:
                            Q_ref = Q_modes_reordered[-1]
                        else:
                            Q_ref = Q_modes_reordered[0]
                    else:
                        Q_ref = None
                    fig = plot_memberships_list(Q_modes_reordered, cmap, names=mode_labels, ind_labels=ind_labels, 
                                                order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label, width_scale=width_scale)
                    fig.savefig(os.path.join(fig_dir,"K{}_modes_{}.png".format(K,suffix)), bbox_inches='tight', dpi=150, transparent=False)
                    plt.close(fig)
            
            if args.plot_type in ["list", "all"]:
                logger.info(f"Plot all modes in a list ({suffix})")
                Q_modes_list = unnest_list(Q_rep_modes_list) if args.use_rep else unnest_list(Q_avg_modes_list)
                mode_names = unnest_list(mode_names_list)
                Q_modes_reordered = [aligned_Qs_allK[mode_name] for mode_name in mode_names]
                mode_labels = [f"{mode_name} ({mode_sizes[mode_name]})" for mode_name in mode_names]
                if args.reorder_ind:
                    if args.reorder_by_max_k:
                        Q_ref = Q_rep_modes_list[-1][0] if args.use_rep else Q_avg_modes_list[-1][0]
                    else:
                        Q_ref = Q_modes_list[0]
                else:
                    Q_ref = None
                fig = plot_memberships_list(Q_modes_reordered, cmap, names=mode_labels, ind_labels=ind_labels, 
                                            order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label, width_scale=width_scale)
                fig.savefig(os.path.join(fig_dir,"all_modes_list_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
                plt.close(fig)

            if args.plot_type in ["major", "all"]:
                logger.info(f"Plot major modes in a list ({suffix})")
                major_mode_names = [mode_names[0] for mode_names in mode_names_list]
                Q_major_modes_reordered = [aligned_Qs_allK[mode_name] for mode_name in major_mode_names]
                major_mode_labels = [f"{mode_name} ({mode_sizes[mode_name]})" for mode_name in major_mode_names]
                if args.reorder_ind:
                    if args.reorder_by_max_k:
                        Q_ref = Q_major_modes_reordered[-1]
                    else:
                        Q_ref = Q_major_modes_reordered[0]
                else:
                    Q_ref = None
                fig = plot_memberships_list(Q_major_modes_reordered, cmap, names=major_mode_labels, ind_labels=ind_labels, 
                                            order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label, width_scale=width_scale)
                fig.savefig(os.path.join(fig_dir,"major_modes_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
                plt.close(fig)      
        else:
            logger.info(f"Plot modes within-K (single K={K_range[0]})")
            i_K = 0
            K = K_range[i_K]
            Q_modes_list = Q_rep_modes_list[i_K] if args.use_rep else Q_avg_modes_list[i_K]
            mode_names = mode_names_list[i_K]
            Q_modes_reordered = reorderQ_within_k(Q_modes_list, mode_names, alignment_acrossK)        
            # Q_modes = cd_res[0]['repQ_modes'] if args.use_rep else cd_res[0]['avgQ_modes']
            suffix = "rep" if args.use_rep else "avg"
            mode_labels = [f"{mode_name} ({mode_sizes[mode_name]})" for mode_name in mode_names]
            if args.reorder_ind:
                if args.reorder_by_max_k:
                    Q_ref = Q_modes_reordered[-1]
                else:
                    Q_ref = Q_modes_reordered[0]
            else:
                Q_ref = None
            fig = plot_memberships_list(Q_modes_reordered, cmap, names=mode_labels, ind_labels=ind_labels, 
                                        order_refQ=Q_ref, order_cls_by_label=args.order_cls_by_label, width_scale=width_scale)
            fig.savefig(os.path.join(fig_dir,"K{}_modes_{}.png".format(K,suffix)), bbox_inches='tight', dpi=150, transparent=False)
            plt.close(fig)

    logger.info(f"Completed".center(50, '-'))
    logger.info(f"".center(50, '='))   

    # zip all files
    logger.info(f"Zipping outputs".center(50, '-'))
    shutil.make_archive(args.output, "zip", args.output)
    logger.info(f"".center(50, '='))  
    tot_toc = time.time()
    logger.info("Total Time: %.3fs", tot_toc-tot_tic)

def parse_args():
    parser = argparse.ArgumentParser(description="Clumppling: a tool for cluster matching and permutation program with integer linear programming")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # required
    required.add_argument("-i", "--input", type=str, required=True, help="Input file path")
    required.add_argument("-o", "--output", type=str, required=True, help="Output file directory")
    required.add_argument("-f", "--format", type=str, required=True, choices=["generalQ", "admixture", "structure", "fastStructure"],
                        help="File format")
    # optional
    optional.add_argument('-v', '--vis', type=str2bool, default=True, required=False, help='Whether to generate figure(s): True (default)/False')
    optional.add_argument('--custom_cmap', type=str, default='', required=False, help='Customized colormap as a comma-separated string of hex codes for colors: if empty (default), using the default colormap, otherwise use the user-specified colormap')
    optional.add_argument("--plot_type", type=str, default="graph", required=False, choices=["graph", "list", "withinK", "major", "all"],
                          help="Type of plot to generate: 'graph' (default), 'list', 'withinK', 'major', 'all'")
    optional.add_argument("--include_cost", type=str2bool, default=True, required=False, help="Whether to include cost values in the graph plot: True (default)/False")
    optional.add_argument("--include_label", type=str2bool, default=True, required=False, help="Whether to include individual labels in the plot: True (default)/False")
    optional.add_argument("--ind_labels", type=str, default="", required=False, 
                        help="A plain text file containing individual labels (one per line) (default: last column from labels in input file, which consists of columns [0, 1, 3] separated by delimiter)")
    optional.add_argument("--regroup_ind", type=str2bool, default=True, required=False, 
                        help="Whether to regroup individuals so that those with the same labels stay together (if labels are available): True (default)/False")
    optional.add_argument("--reorder_ind", type=str2bool, default=True, required=False, 
                        help="Whether to reorder individuals within each label group in the plot (if labels are available): True (default)/False")
    optional.add_argument("--reorder_by_max_k", type=str2bool, default=True, required=False, 
                        help="Whether to reorder individuals based on the major mode with largest K: True (default)/False (based on the major mode with smallest K)")
    optional.add_argument("--order_cls_by_label", type=str2bool, default=True, required=False, 
                        help="Whether to reorder clusters based on total memberships within each label group in the plot: True (default)/False (by overall total memberships)")

    optional.add_argument("--extension", type=str, default="", required=False, help="Extension of input files")
    optional.add_argument("--skip_rows", type=int, default=0, required=False, help="Skip top rows in input files")
    optional.add_argument("--remove_missing", type=str2bool, default=True, required=False, help="Remove individuals with missing data: True (default)/False")

    optional.add_argument("--cd_method", type=str, default="louvain", required=False,
                          choices=["louvain", "leiden", "infomap", "markov_clustering", "label_propagation", "walktrap", "custom"],
                          help="Community detection method to use (default: louvain)")
    optional.add_argument("--cd_res", type=float, default=1.0, required=False,
                          help="Resolution parameter for the default Louvain community detection (default: 1.0)")
    optional.add_argument("--test_comm", type=str2bool, default=True, required=False,
                          help="Whether to test community structure (default: True)")
    optional.add_argument("--comm_min", type=float, default=1e-4, required=False,
                          help="Minimum threshold for cost matrix (default: 1e-4)")
    optional.add_argument("--comm_max", type=float, default=1e-2, required=False,
                          help="Maximum threshold for cost matrix (default: 1e-2)")
    optional.add_argument("--merge", type=str2bool, default=True, required=False,
                          help="Whether to merge two clusters when aligning K+1 to K (default: True)")
    optional.add_argument("--use_rep", type=str2bool, default=True, required=False, help="Use representative modes (alternative: average): True (default)/False")
    optional.add_argument("--use_best_pair", type=str2bool, default=True, required=False, help="Use best pair as anchor for across-K alignment (alternative: major): True (default)/False")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(args.output, "clumppling.log"))
    disp_params(args, title="CLUMPPLING")
    main(args)