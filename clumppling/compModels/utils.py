import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from clumppling.plot import plot_membership

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
    mode_stats['SizeFrac'] = mode_stats.apply(lambda row: f"{row['Size']}/{size_K[int(row['Mode'].split("M")[0].strip("K"))]}", axis=1)

    return mode_stats


def plot_multi_model_graph(K_range: list[int], models: list[str],
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
    n_K_in_model = np.zeros((len(K_range),len(models)))
    for i_K, K in enumerate(K_range):
        mode_names = Q_names_list[i_K]
        for i_model, model in enumerate(models):
            n_K_in_model[i_K, i_model] = len([name for name in mode_names if name.startswith(model)])
    n_row = np.sum(n_K_in_model.max(axis=1)).astype(int) 
    n_col = len(models)

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
                ax = fig.add_subplot(gs[np.sum(n_K_in_model[:i_K,i_model]).astype(int)+i_mode, i_model])
                plot_membership(Q_mats[i_mode], cmap, ax=ax, title=mode_name_labels_model_K[i_mode])
    # label each K
    if label_K:
        K_row_indices = np.insert(np.cumsum(n_K_in_model.max(axis=1)).astype(int),0,0)[:-1]
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
            fig.text(x_center, 0.92, f"{models[c].title()}", 
                    va="top", ha="center", fontsize=16, fontweight="bold")
    return fig