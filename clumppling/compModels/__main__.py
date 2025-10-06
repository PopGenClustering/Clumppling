import os
import numpy as np
import shutil
import argparse
import matplotlib.pyplot as plt

from clumppling.log_config import setup_logger
from clumppling.utils import disp_params, disp_msg, str2bool, parse_strings, unnest_list, write_reordered_across_k
from clumppling.core import align_across_k, write_alignment_across_k, reorderQ_across_k
from clumppling.plot import plot_colorbar, load_default_cmap, parse_custom_cmap, plot_alignment_list, plot_alignment_graph
from clumppling.alignAcrossK import load_any_qfiles, extract_K_range_from_Qs, separate_Qs_by_K
from . import plot_multi_model_list, load_mode_stats, plot_multi_model_graph_sbs, plot_multi_model_graph_il

import logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="clumppling.compModels")

    parser.add_argument("--models", nargs="+", type=str, required=True,
                        help="List of model names.")
    parser.add_argument("--qfilelists", nargs="+", type=str, required=True,
                        help="List of files containing Q file names from each model.")
    parser.add_argument("--qnamelists", nargs="+", type=str, default="",
                        help="List of files containing replicate names from each model.")
    parser.add_argument("--mode_stats_files", nargs="+", type=str, default="",
                        help="List of files containing mode statistics from each model.")
    parser.add_argument("--ind_labels", type=str, default="", required=False, 
                        help="A plain text file containing individual labels (one per line)")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Output file directory")
    parser.add_argument('-v', '--vis', type=str2bool, default=True, required=False, help='Whether to generate figure(s): True (default)/False')
    parser.add_argument('--custom_cmap', type=str, default='', required=False, help='A plain text file containing customized colors (one per line; in hex code): if empty (default), using the default colormap, otherwise use the user-specified colormap')
    parser.add_argument('--include_sim_in_label', type=str2bool, default=True, required=False, help='Whether to include (original) alignment similarity in mode labels (if provided): True (default)/False')
    parser.add_argument("--fig_format", type=str, default="tiff", required=False, choices=["png", "jpg", "jpeg", "tif", "tiff", "svg", "pdf", "eps", "ps", "bmp", "gif"], help="Figure format for output files (default: tiff)")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logger(os.path.join(args.output, "compModels.log"))
    disp_params(args, title="Comparison of Models")

    if os.path.exists(args.output) and os.listdir(args.output):
        shutil.rmtree(args.output)
        logger.info(f"Output directory '{args.output}' already exists and is not empty. Removed existing directory.")
    os.makedirs(args.output, exist_ok=True)
    logger.info(f"Created output directory '{args.output}'.")

    assert len(args.qfilelists)==len(args.models) >=2, "At least two models are required for comparison."
    disp_msg(f"Loading replicate names and matrices from each model")
    Q_files_list = [parse_strings(f) for f in args.qfilelists]
    if args.qnamelists is None:
        Q_names_list = [[os.path.splitext(os.path.basename(p))[0] for p in Q_files] for Q_files in Q_files_list]
    else:
        Q_names_list = [parse_strings(f) for f in args.qnamelists]
    Q_mat_list = [load_any_qfiles(Q_files) for Q_files in Q_files_list]

    if args.ind_labels!="":
        ind_labels = parse_strings(args.ind_labels, remove_dup=False)
        disp_msg(f"Loaded individual labels for {len(ind_labels)} individuals from {args.ind_labels}")
    
    for i_model, model in enumerate(args.models):
        if "_" in model or "-" in model:
            raise ValueError(f"Model name '{model}' should not contain underscore or dash.")
    # new names for replicates
    Q_new_names_list = list()
    for i_model, model in enumerate(args.models):
        Q_new_names_list.append([f"{model}_{name}" for name in Q_names_list[i_model]])

    Q_list = unnest_list(Q_mat_list)
    Q_names = unnest_list(Q_new_names_list)
    K_range = extract_K_range_from_Qs(Q_list)
    Q_modes_list, mode_names_list = separate_Qs_by_K(K_range, Q_list, Q_names)
    alignment_acrossK, cost_acrossK, best_acrossK_out, major_acrossK_out = align_across_k(K_range, Q_modes_list, mode_names_list, merge=True)
    disp_msg(f"Saving alignment results to {args.output}")
    write_alignment_across_k(alignment_acrossK, cost_acrossK, os.path.join(args.output,"alignment_across_all.txt"))
    best_acrossK_out.to_csv(os.path.join(args.output,"best_pairs_acrossK.txt"), index=False)
    
    for i_model, model in enumerate(args.models):
        Q_new_names = Q_new_names_list[i_model]
        np.savetxt(os.path.join(args.output, f"{model}_replicates.txt"),
                   Q_new_names, fmt="%s")
    
    # save aligned replicates
    disp_msg(f"Reordering membership matrices according to alignment")
    modeQ_dir = os.path.join(args.output,"aligned")
    anchor_pairs = best_acrossK_out["Best Pair"].tolist()
    aligned_Qs_allK, all_modes_alignment = reorderQ_across_k(K_range, Q_modes_list, mode_names_list, 
                                                            alignment_acrossK, anchor_pairs)
    write_reordered_across_k(aligned_Qs_allK, all_modes_alignment, output_dir=modeQ_dir, suffix="")
    
    disp_msg(f"Plotting alignment results")
    if args.vis:
        K_max = max(K_range)
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
        logger.info(f"Plot alignment patterns")
        Q_names_reordered = unnest_list(mode_names_list)
        Q_aligned_list = load_any_qfiles([os.path.join(modeQ_dir,f"{Q_name}.Q") for Q_name in Q_names_reordered])
        mode_K = [Q.shape[1] for Q in Q_aligned_list]
        fig = plot_alignment_list(mode_K, Q_names_reordered, cmap, alignment_acrossK, all_modes_alignment, marker_size=250)
        fig.axes[0].set_ylabel("Replicates")
        fig.savefig(os.path.join(fig_dir,"alignment_pattern_list.{}".format(args.fig_format)), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig)

        # plot alignment pattern graph 
        y_aspect = 3 if K_max<5 else 3+K_max/5
        wspace_padding = 1.3 if K_max<8 else 1.15
        fig = plot_alignment_graph(K_range, names_list=mode_names_list, cmap=cmap, 
                                   alignment_acrossK=alignment_acrossK, all_modes_alignment=all_modes_alignment,
                                   anchor_pairs=anchor_pairs,
                                   alt_color=False, ls_alt=['-', '--'],  
                                   y_aspect=y_aspect, wspace_padding=wspace_padding, separate_labels=True) #color_alt=['#6A8A9F','#B49F63', '#789B8C','#AD839A','#8D8D8D'], 
        fig.savefig(os.path.join(fig_dir,"alignment_pattern_graph.{}".format(args.fig_format)), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig)

        # load mode statistics if provided
        if args.mode_stats_files and len(args.mode_stats_files)==len(args.models):
            disp_msg(f"Loading mode statistics from each model")
            mode_stats_list = [load_mode_stats(f) for f in args.mode_stats_files]
            mode_sizes = dict()
            mode_performances = dict()
            for i_model, model in enumerate(args.models):
                mode_stats = mode_stats_list[i_model]
                if 'Mode' not in mode_stats.columns or 'Size' not in mode_stats.columns or 'Performance' not in mode_stats.columns:
                    logger.warning(f"Mode statistics file '{args.mode_stats_files[i_model]}' does not contain required columns 'Mode', 'Size', and 'Performance'. Skipping mode size labeling for model '{model}'.")
                    continue
                for _, row in mode_stats.iterrows():
                    mode_name = f"{model}_{row['Mode']}"
                    if 'SizeFrac' in mode_stats.columns:
                        mode_sizes[mode_name] = row['SizeFrac']
                    else:
                        mode_sizes[mode_name] = row['Size']
                    mode_performances[mode_name] = row['Performance']

        if args.include_sim_in_label and 'Performance' in mode_stats.columns:
            mode_labels_list = [["{} ({}): sim {:.3f}".format(mode_name.replace('_', ' '), mode_sizes[mode_name], mode_performances[mode_name]) for mode_name in mode_names] for mode_names in mode_names_list]
        else:
            mode_labels_list = [[f"{mode_name.replace('_', ' ')} ({mode_sizes[mode_name]})" for mode_name in mode_names] for mode_names in mode_names_list]        
        
        logger.info(f"Plot structure plots of aligned modes from all models")
        # plot graph of structure plots interleaved
        fig = plot_multi_model_graph_il(K_range, args.models, 
                                     mode_names_list, Q_names_reordered, 
                                     Q_aligned_list, cmap, ind_labels=ind_labels,
                                     Q_names_label_list=mode_labels_list,
                                     label_K=True, label_model=True)
        fig.savefig(os.path.join(fig_dir,"comparison_aligned_models_graph_il.{}".format(args.fig_format)), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig) 

        # plot graph of structure plots side by side
        fig = plot_multi_model_graph_sbs(K_range, args.models, 
                                     mode_names_list, Q_names_reordered, 
                                     Q_aligned_list, cmap, 
                                     Q_names_label_list=mode_labels_list,
                                     label_K=True, label_model=True)
        fig.savefig(os.path.join(fig_dir,"comparison_aligned_models_graph_sbs.{}".format(args.fig_format)), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig) 
        
        # plot list of structure plots side by side
        fig = plot_multi_model_list(K_range, args.models, 
                                     mode_names_list, Q_names_reordered, 
                                     Q_aligned_list, cmap, 
                                     Q_names_label_list=mode_labels_list,
                                     label_K=True, label_model=True)
        fig.savefig(os.path.join(fig_dir,"comparison_aligned_models_list.png"), bbox_inches='tight', dpi=150, transparent=False)
        plt.close(fig)  

    # zip all files
    logger.info(f"Zipping outputs".center(50, '-'))
    shutil.make_archive(args.output, "zip", args.output)

    # complete
    logger.info(f"Completed".center(50, '-'))
    logger.info(f"".center(50, '=')) 