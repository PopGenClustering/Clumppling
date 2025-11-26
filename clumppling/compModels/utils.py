import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from clumppling.plot import plot_membership
from clumppling.utils import get_uniq_lb_sep
import logging
logger = logging.getLogger(__name__)

def load_mode_stats(mode_stats_file: str) -> pd.DataFrame:
    """
    Load mode statistics from a file.
    Args:
        mode_stats_file (str): Path to the mode statistics file.

    Returns:
        pd.DataFrame: A DataFrame containing mode statistics.
    """
    mode_stats = pd.read_csv(mode_stats_file, sep=",", header=0)
    # set size to s/total over each K
    size_K = dict()
    for _, row in mode_stats.iterrows():
        K = int(row['Mode'].split("M")[0].strip("K"))  # extract K from mode name
        if K in size_K:
            size_K[K] += row['Size']
        else:
            size_K[K] = row['Size']
    mode_stats['SizeFrac'] = mode_stats.apply(lambda row: f"{row['Size']}/{size_K[int(row['Mode'].split('M')[0].strip('K'))]}", axis=1)

    return mode_stats


def plot_multi_model_list(K_range: list[int], models: list[str],
                           Q_names_list: list[list[str]],
                           Q_names_reordered: list[str],
                           Q_aligned_list: list[np.ndarray],
                           cmap: list[tuple[float, float, float]],
                           Q_names_label_list: list[list[str]] = [],
                           label_K: bool = True,
                           label_model: bool = True):
    """
    Plot membership matrices for multiple models across different K values.
    Args:
        K_range (list[int]): List of K values.
        models (list[str]): List of model names.
        Q_names_list (list[list[str]]): List of Q names for each K.
        Q_names_reordered (list[str]): Reordered list of Q names.
        Q_aligned_list (list[np.ndarray]): List of aligned Q matrices.
        cmap (list[tuple[float, float, float]]): Colormap for plotting.
        label_K (bool): Whether to label each K value. Default is True.
        label_model (bool): Whether to label each model. Default is True.
    """
    if len(Q_names_label_list) == 0:
        Q_names_label_list = Q_names_list
    # get number of replicates for each K for each model
    n_row = 0
    n_K_in_model = np.zeros((len(K_range),len(models)), dtype=int)
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        for i_model, model in enumerate(models):
            n_K_in_model[i_K, i_model] = len([name for name in mode_names if name.startswith(model)])
    n_row = np.sum(n_K_in_model.max(axis=1))
    n_col = len(models)
    idx_K_row_start = np.insert(np.cumsum(n_K_in_model.max(axis=1)),0,0)[:-1]

    fig = plt.figure(figsize=(7*n_col,2*n_row), dpi=150)
    gs = fig.add_gridspec(n_row,n_col, wspace=0.2, hspace=0.5)
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        mode_name_labels = Q_names_label_list[i_K]
        for i_model, model in enumerate(models):
            Q_names_model_K = [name for name in mode_names if name.startswith(model)]
            mode_name_labels_model_K = [mode_name_labels[i_name] for i_name,name in enumerate(mode_names) if name.startswith(model)]
            # get corresponding Q matrices
            Q_mats = [Q_aligned_list[Q_names_reordered.index(name)] for name in Q_names_model_K]
            for i_mode in range(len(Q_mats)):
                row_idx = idx_K_row_start[i_K] + i_mode
                ax = fig.add_subplot(gs[row_idx, i_model])
                plot_membership(Q_mats[i_mode], cmap, ax=ax, title=mode_name_labels_model_K[i_mode])
    # label each K
    if label_K:
        K_row_indices = np.insert(np.cumsum(n_K_in_model.max(axis=1)),0,0)[:-1]
        for row in range(gs.nrows):
            if row in K_row_indices:
                # Get the GridSpec cell rectangle in figure coordinates
                bb = gs[row, 0].get_position(fig)
                y_center = (bb.y0 + bb.y1) / 2  # vertical midpoint
                fig.text(0.07, y_center, f"K={K_range[np.where(K_row_indices==row)[0][0]]}", fontsize=16, weight='bold', rotation=0, va="center", ha="left")
    # label each column (model)
    if label_model:
        for c in range(gs.ncols):
            bb = gs[0, c].get_position(fig)  # position of first row in that column
            x_center = (bb.x0 + bb.x1) / 2
            fig.text(x_center, 0.92, f"{models[c]}", 
                    va="top", ha="center", fontsize=16, fontweight="bold")
    return fig


def plot_multi_model_graph_sbs(K_range: list[int], models: list[str],
                           Q_names_list: list[list[str]],
                           Q_names_reordered: list[str],
                           Q_aligned_list: list[np.ndarray],
                           cmap: list[tuple[float, float, float]],
                           Q_names_label_list: list[list[str]] = [],
                           label_K: bool = True,
                           label_model: bool = True):
    """
    Plot membership matrices for multiple models across different K values.
    Args:
        K_range (list[int]): List of K values.
        models (list[str]): List of model names.
        Q_names_list (list[list[str]]): List of Q names for each K.
        Q_names_reordered (list[str]): Reordered list of Q names.
        Q_aligned_list (list[np.ndarray]): List of aligned Q matrices.
        cmap (list[tuple[float, float, float]]): Colormap for plotting.
        label_K (bool): Whether to label each K value. Default is True.
        label_model (bool): Whether to label each model. Default is True.
    """
    if len(Q_names_label_list) == 0:
        Q_names_label_list = Q_names_list
    # get number of replicates for each K for each model
    n_K_in_model = np.zeros((len(K_range),len(models)), dtype=int)
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        for i_model, model in enumerate(models):
            n_K_in_model[i_K, i_model] = len([name for name in mode_names if name.startswith(model)])
    n_row = len(K_range)
    col_in_models = n_K_in_model.max(axis=0)
    n_col = np.sum(col_in_models)

    fig = plt.figure(figsize=(7*n_col,2*n_row), dpi=150)
    gs = fig.add_gridspec(n_row,n_col, wspace=0.2, hspace=0.5)
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        mode_name_labels = Q_names_label_list[i_K]

        for i_model, model in enumerate(models):
            Q_names_model_K = [name for name in mode_names if name.startswith(model)]
            mode_name_labels_model_K = [mode_name_labels[i_name] for i_name,name in enumerate(mode_names) if name.startswith(model)]
            # get corresponding Q matrices
            Q_mats = [Q_aligned_list[Q_names_reordered.index(name)] for name in Q_names_model_K]
            for i_mode in range(len(Q_mats)):
                ax = fig.add_subplot(gs[i_K, sum(col_in_models[:i_model])+i_mode])
                plot_membership(Q_mats[i_mode], cmap, ax=ax, title=mode_name_labels_model_K[i_mode])
    # label each K
    if label_K:
        for row in range(gs.nrows):
            # Get the GridSpec cell rectangle in figure coordinates
            bb = gs[row, 0].get_position(fig)
            y_center = (bb.y0 + bb.y1) / 2  # vertical midpoint
            fig.text(0.11, y_center, f"K={K_range[row]}", 
                     fontsize=17, weight='bold', rotation=0, va="center", ha="left")
    # label each column (model)
    if label_model:
        for i_c,c in enumerate(np.insert(np.cumsum(col_in_models), 0, 0)[:-1]):
            bb = gs[0, c].get_position(fig)  # position of first row in that column
            x_left = bb.x0
            fig.text(x_left, 0.91, f"{models[i_c]}", 
                    va="bottom", ha="left", fontsize=17, fontweight="bold")
    return fig


def plot_multi_model_graph_il(K_range: list[int], models: list[str],
                           Q_names_list: list[list[str]],
                           Q_names_reordered: list[str],
                           Q_aligned_list: list[np.ndarray],
                           cmap: list[tuple[float, float, float]],
                           ind_labels: list[str] = [],
                           Q_names_label_list: list[list[str]] = [],
                           label_K: bool = True,
                           label_model: bool = True, 
                           bg_colors: list = []):
    """
    Plot membership matrices for multiple models across different K values.
    Args:
        K_range (list[int]): List of K values.
        models (list[str]): List of model names.
        Q_names_list (list[list[str]]): List of Q names for each K.
        Q_names_reordered (list[str]): Reordered list of Q names.
        Q_aligned_list (list[np.ndarray]): List of aligned Q matrices.
        cmap (list[tuple[float, float, float]]): Colormap for plotting.
        label_K (bool): Whether to label each K value. Default is True.
        label_model (bool): Whether to label each model. Default is True.
    """
    if len(Q_names_label_list) == 0:
        Q_names_label_list = Q_names_list
    if len(ind_labels)>0:
        uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx = get_uniq_lb_sep(ind_labels)
        
    # get number of replicates for each K for each model
    n_K_in_model = np.zeros((len(K_range),len(models)), dtype=int)
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        for i_model, model in enumerate(models):
            n_K_in_model[i_K, i_model] = len([name for name in mode_names if name.startswith(model)])
    n_models = len(models)
    n_row = n_models*len(K_range)
    n_col = np.max(n_K_in_model)

    if len(bg_colors) < n_models:
        gray_cmap = cm.get_cmap('gray')
        bg_colors = [gray_cmap(point) for point in np.linspace(0.85, 1, n_models)]
    elif len(bg_colors) >= n_models:
        bg_colors = bg_colors[:n_models]

    fig = plt.figure(figsize=(7*n_col,2*n_row), dpi=150)
    gs = fig.add_gridspec(n_row,n_col, wspace=0.2, hspace=0.5)

    for row in range(n_row):
        ss = gs[row, :]                # row across all columns
        pos = ss.get_position(fig)
        color = bg_colors[row % n_models]

        rect = Rectangle(
            (pos.x0-0.01, pos.y0-0.01),
            pos.width + 0.02,
            pos.height + 0.02,
            transform=fig.transFigure,
            color=color,
            zorder=0,
            clip_on=False, 
        )
        fig.add_artist(rect)

    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        mode_name_labels = Q_names_label_list[i_K]

        for i_model, model in enumerate(models):
            Q_names_model_K = [name for name in mode_names if name.startswith(model)]
            mode_name_labels_model_K = [mode_name_labels[i_name] for i_name,name in enumerate(mode_names) if name.startswith(model)]
            # get corresponding Q matrices
            Q_mats = [Q_aligned_list[Q_names_reordered.index(name)] for name in Q_names_model_K]
            for i_mode in range(len(Q_mats)):
                ax = fig.add_subplot(gs[i_K*n_models + i_model, i_mode])

                plot_membership(Q_mats[i_mode], cmap, ax=ax, title=mode_name_labels_model_K[i_mode])
                if len(ind_labels)>0:
                    for v in uniq_lbs_sep_idx:
                        ax.axvline(v, ymin=-0.2, ymax=1, color='k', ls='--', lw=0.5, clip_on=False)
                if n_K_in_model[i_K,i_model]>0 and i_K==len(K_range)-1 and i_model==len(models)-1: #or (n_K_in_model[i_K,i_model]>0 and n_K_in_model[i_K+1,i_model]==0)
                    ax.set_xticks(uniq_lbs_indices)
                    ax.tick_params(axis='x', which='both', length=0, pad=2)
                    if len(uniq_lbs)>=10:
                        rot = 90
                    else:
                        rot = 45 if np.any([len(lb)>5 for lb in uniq_lbs]) else 0
                    fs = 14 if len(uniq_lbs)<10 else 10
                    ax.set_xticklabels(uniq_lbs, rotation=rot, ha='center', fontsize=fs)
                    # ax.tick_params(axis="x", pad=2)
    # label each K
    if label_K:
        for row in range(gs.nrows):
            i_model = row % n_models
            i_K = row // n_models
            # Get the GridSpec cell rectangle in figure coordinates
            bb = gs[row, 0].get_position(fig)
            y_center = (bb.y0 + bb.y1) / 2  # vertical midpoint
            fig.text(bb.x0-0.015, y_center, f"{models[i_model]} K={K_range[i_K]}", 
                        fontsize=17, weight='bold', rotation=0, va="center", ha="right")
    # label each column (model)
    if label_model:
        for col in range(gs.ncols):
            bb = gs[-1, col].get_position(fig)  # position of first row in that column
            fig.text((bb.x0+bb.x1)/2, bb.y0-0.015, f"M{col+1}", 
                    va="top", ha="center", fontsize=17, fontweight="bold")
    return fig