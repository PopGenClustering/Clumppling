import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
import matplotlib.colors as mcolors
import colorsys
import importlib.resources
from matplotlib import cm
from matplotlib.patches import ConnectionPatch, Rectangle
from typing import Optional, cast
from clumppling.utils import get_uniq_lb_sep, reorder_ind_within_group

import logging
logger = logging.getLogger(__name__)

def load_default_cmap(K:int):
    """
    Load the default color map from the resource file.
    
    Returns:
        A list of colors (in hex code) loaded from the resource file.
    """
    resource_file = importlib.resources.files('clumppling.files').joinpath('default_colors.txt')
    with resource_file.open('r') as f:
        colors = f.readlines()
    colors = [color.strip() for color in colors if color.strip()]  # Clean up empty lines
    cmap = mcolors.ListedColormap(colors[:K])
    rgb_list = [mcolors.to_rgb(c) for c in np.array(cmap.colors)]
    return rgb_list

def parse_custom_cmap(colors: list, K: int) -> list[tuple[float, float, float]]:
    """
    Parse a custom color map string into a list of colors.
    
    Args:
        cmap_str: A string of comma-separated color codes (e.g., "#FF0000,#00FF00,#0000FF").
        K: The number of colors to return.
        
    Returns:
        A list of colors (in hex code) parsed from the input string.
    """
    # colors = [color.strip() for color in cmap_str.split(',') if color.strip()]
    if len(colors) < K:
        raise ValueError(f"WARNING: Custom color map has only {len(colors)} colors, but K={K} is requested. Recycling colors.")
    cmap = mcolors.ListedColormap(colors, N=K)
    rgb_list = [mcolors.to_rgb(c) for c in np.array(cmap.colors)]
    return rgb_list


def plot_colorbar(rgb_colors: list[tuple[float, float, float]], K_max: int, fig_dir:str, fig_format: str="tiff"):
    assert K_max == len(rgb_colors), f"Number of colors in colormap ({len(rgb_colors)}) is K_max ({K_max})"

    fig, ax = plt.subplots(figsize=(K_max, 1),facecolor='white',subplot_kw=dict(yticks=[])) 
    ax.set_xticks(ticks=np.arange(0.5,K_max+0.5,1),labels=["Cls.{}".format(k) for k in np.arange(1,K_max+1)])
    img = np.array([rgb_colors])
    ax.imshow(img, aspect='auto', extent=(0, K_max, 0, 1))
    ax.get_yaxis().set_visible(False)

    fig.savefig(os.path.join(fig_dir,"colorbar.{}".format(fig_format)), bbox_inches='tight',dpi=300, transparent=False)
    plt.close(fig)


def plot_membership(Q:np.ndarray, cmap: list[tuple[float, float, float]], 
                    ax: Optional[maxes.Axes] = None, ylab: str="", title: str="",
                    fontsize: float=14) -> maxes.Axes:
    if ax is None:
        _, ax = plt.subplots()
        ax = cast(maxes.Axes, ax)

    N = Q.shape[0]
    K = Q.shape[1]
    Q_aug = np.hstack((np.zeros((N,1)),Q))
    for k in range(K):
        ax.bar(range(N), Q_aug[:,(k+1)], bottom=np.sum(Q_aug[:,0:(k+1)],axis=1), 
               width=1.0, edgecolor='w', linewidth=0, facecolor=cmap[k])
    ax.set_xlim(0,N)
    ax.set_xticks([])
    ax.set_ylim(0,1)
    ax.set_yticks([])
    if ylab:
        ax.set_ylabel("\n".join(ylab.split()), rotation=0, fontsize=fontsize, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    if title:
        ax.set_title(title, fontsize=fontsize, loc="left", pad=5)
    return ax


def plot_membership_reordered(Q:np.ndarray, lbs: list, 
                              mbsp_sort_indices: dict,
                              cmap: list[tuple[float, float, float]], 
                              ax: Optional[maxes.Axes] = None, ylab: str="", title: str="", 
                              fontsize: float=14) -> maxes.Axes:
    if ax is None:
        _, ax = plt.subplots()
        ax = cast(maxes.Axes, ax)
    N = Q.shape[0]
    K = Q.shape[1]

    for lb in np.unique(lbs):
        lb_indices = [i for i, val in enumerate(lbs) if val == lb]
        lb_P = Q[lb_indices,:]
        mbsp_sortidx = mbsp_sort_indices[lb]
        lb_P_sorted = lb_P[:,mbsp_sortidx]
        lb_P_aug = np.hstack((np.zeros((lb_P_sorted.shape[0],1)),lb_P_sorted))
        for k in range(K):
            ax.bar(lb_indices, lb_P_aug[:,(k+1)], bottom=np.sum(lb_P_aug[:,0:(k+1)],axis=1), 
                    width=1.0, edgecolor='w', linewidth=0, facecolor=cmap[mbsp_sortidx[k]])
    ax.set_xlim(0,N)
    ax.set_xticks([])
    ax.set_ylim(0,1)
    ax.set_yticks([])
    if ylab:
        ax.set_ylabel("\n".join(ylab.split()), rotation=0, fontsize=fontsize, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    if title:
        ax.set_title(title, fontsize=fontsize, loc="left", pad=5)
    return ax


def plot_memberships_list(Q_list: list[np.ndarray], cmap: list[tuple[float, float, float]], 
                          names: list[str], ind_labels: list=[], right_labels: Optional[list[str]]=[],
                          order_refQ: Optional[np.ndarray]=None, order_cls_by_label: bool=True,
                          width_scale: float=1.0):
    # get separation points in the labels
    if len(ind_labels)>0:
        uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx = get_uniq_lb_sep(ind_labels)
        if order_refQ is not None:
            ref_ind_sorted_indices, ref_mbsp_sort_indices = reorder_ind_within_group(order_refQ, ind_labels)
            # print(f"Reordering individuals within each label group based on the reference Q matrix with shape {order_refQ.shape}.")

    fig, axes = plt.subplots(len(Q_list),1,figsize=(15*width_scale,1.5*len(Q_list)),facecolor='white')

    for i in range(len(Q_list)):
        Q = Q_list[i]
        lb = names[i]
        right_lb = right_labels[i] if right_labels and i<len(right_labels) else ""
        ax = axes[i] if len(Q_list)>1 else axes

        if order_refQ is None or len(ind_labels)==0:
            plot_membership(Q, cmap, ax=ax, ylab=lb)
        else:
            if order_cls_by_label:
                _, mbsp_sort_indices = reorder_ind_within_group(Q, ind_labels)
                reordered_ind_labels = [ind_labels[i] for i in ref_ind_sorted_indices]
                plot_membership_reordered(Q[ref_ind_sorted_indices, :], reordered_ind_labels, mbsp_sort_indices, cmap, ax=ax, ylab=lb)
            else:
                Q_reordered = Q[ref_ind_sorted_indices, :]
                mbsp_sum = Q_reordered.sum(axis=0)
                mbsp_sortidx = np.argsort(-mbsp_sum) # from largest to smallest
                Q_reordered = Q_reordered[:,mbsp_sortidx]
                plot_membership(Q_reordered, [cmap[i] for i in mbsp_sortidx], ax=ax, ylab=lb)        
        ax.set_title(right_lb, fontsize=14, loc="right")

        if len(ind_labels)>0:
            if len(ind_labels)!=Q.shape[0]:
                logger.warning(f"Number of individual labels ({len(ind_labels)}) does not match number of rows in Q matrix ({Q.shape[0]}). ")
            for v in uniq_lbs_sep_idx:
                ax.axvline(v, ymin=-0.2, ymax=1, color='k', ls='--', lw=0.5, clip_on=False)
            if i==len(Q_list)-1:
                ax.set_xticks(uniq_lbs_indices)
                ax.tick_params(axis='x', which='both', length=0)
                rot = 45 if np.any([len(lb)>5 for lb in uniq_lbs]) else 0
                fs = 14 if len(uniq_lbs)<10 else 10
                ax.set_xticklabels(uniq_lbs, rotation=rot, ha='center', fontsize=fs)
    fig.subplots_adjust(hspace=0.3)
    return fig

def adjust_lightness(color, factor=1.0):
    """
    Adjust color lightness by multiplying the lightness channel in HLS.
    factor > 1 → lighter, factor < 1 → darker
    """
    r, g, b = mcolors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0, min(1, l * factor))
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return (r, g, b)


def create_single_cmap(color: str, name: str='custom_cmap', light: float=1.5, dark: float=0.6) -> mcolors.LinearSegmentedColormap:
    """Create a colormap that transitions from white to the specified color and then to black.
    Args:
        color (str): The base color for the colormap.
        name (str): The name of the colormap.
    """
    rgb = mcolors.to_rgb(color)
    # Define gradient: white → base color → black
    return mcolors.LinearSegmentedColormap.from_list(
        name, 
        [adjust_lightness(color, light), mcolors.to_rgb(color), adjust_lightness(color, dark)]
    )


def plot_graph(K_range: list[int], Q_list_list: list[list[np.ndarray]], cmap: list[tuple[float, float, float]], 
               names_list: list[list[str]], labels_list: Optional[list[list[str]]]=None, right_labels_list: Optional[list[list[str]]]=None,
               cost_acrossK: Optional[dict]=None,
               ind_labels: list=[], fontsize: float=14, 
               alt_color: bool=True, line_cmap=None,
               color_alt: list[str]=['rosybrown', 'steelblue','goldenrod', 'darkseagreen'], 
               order_refQ: Optional[np.ndarray]=None, order_cls_by_label: bool=True, 
               width_scale: float=1.0, height_scale: float=1.0):
    mode_numbers = [len(names_list[i_K]) for i_K, K in enumerate(K_range)]
    n_row = len(K_range)
    n_col = max(mode_numbers)
    N = Q_list_list[0][0].shape[0]

    # get separation points in the labels
    if len(ind_labels)>0:
        uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx = get_uniq_lb_sep(ind_labels)
        if order_refQ is not None:
            ref_ind_sorted_indices, ref_mbsp_sort_indices = reorder_ind_within_group(order_refQ, ind_labels)

    if labels_list is None:
        labels_list = names_list

    fig = plt.figure(figsize=(7*n_col*width_scale,2.5*n_row*height_scale), dpi=150)
    gs = fig.add_gridspec(n_row,n_col, wspace=0.2, hspace=1)

    axes_handles = dict()
    for i_K, K in enumerate(K_range):
        for i_mode in range(mode_numbers[i_K]):
            ax = fig.add_subplot(gs[i_K, i_mode])
            axes_handles[(K, i_mode)] = ax
            Q_label = labels_list[i_K][i_mode] 
            right_label = right_labels_list[i_K][i_mode] if right_labels_list is not None else ""
            Q = Q_list_list[i_K][i_mode]
            if order_refQ is None or len(ind_labels)==0:
                plot_membership(Q, cmap, ax=ax, title=Q_label, fontsize=fontsize)
            else:
                if order_cls_by_label:
                    _, mbsp_sort_indices = reorder_ind_within_group(Q, ind_labels)
                    reordered_ind_labels = [ind_labels[i] for i in ref_ind_sorted_indices]
                    plot_membership_reordered(Q[ref_ind_sorted_indices, :], reordered_ind_labels, mbsp_sort_indices, cmap, ax=ax, title=Q_label, fontsize=fontsize)
                else:
                    Q_reordered = Q[ref_ind_sorted_indices, :]
                    mbsp_sum = Q_reordered.sum(axis=0)
                    mbsp_sortidx = np.argsort(-mbsp_sum) # from largest to smallest
                    Q_reordered = Q_reordered[:,mbsp_sortidx]
                    plot_membership(Q_reordered, [cmap[i] for i in mbsp_sortidx], ax=ax, title=Q_label, fontsize=fontsize) 
            ax.set_title(right_label, fontsize=14, loc="right")

            ax.set_zorder(9)
            if len(ind_labels)>0:
                if len(ind_labels)!=N:
                    logger.warning(f"Number of individual labels ({len(ind_labels)}) does not match number of rows in Q matrix ({N}). ")
                for v in uniq_lbs_sep_idx:
                    ax.axvline(v, ymin=-0.2, ymax=1, color='k', ls='--', lw=0.5, clip_on=False)
                if i_K==len(K_range)-1:
                    ax.set_xticks(uniq_lbs_indices)
                    ax.tick_params(axis='x', which='both', length=0)
                    rot = 45 if np.any([len(lb)>5 for lb in uniq_lbs]) else 0
                    if len(uniq_lbs)>10:
                        rot = 90 
                    fs = 14 if len(uniq_lbs)<10 else 10
                    ax.set_xticklabels(uniq_lbs, rotation=rot, ha='center', fontsize=fs)

    for i_K, K in enumerate(K_range):
        ax = axes_handles[(K, 0)]
        ax.set_ylabel(f"K={K}", fontsize=18, rotation=0, labelpad=10, va="center", ha="right", weight='bold')

    # plot cost if applicable
    def get_frac_coord(coord1,coord2,frac: float=0.5):
        return (coord2[0]-coord1[0])*frac+coord1[0], (coord2[1]-coord1[1])*frac+coord1[1] 
    
    if cost_acrossK is not None:
        
        logger.info(f"Including graph edges for cost across K")
        textbox_props = dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='none', pad=0.1)
        
        # get colormap for connection lines
        if not alt_color:
            if line_cmap is None:
                logger.info(f"Using default Greys colormap for cost lines")
                line_cmap = cm.get_cmap("Greys")
        else:
            cmaps_alt = [create_single_cmap(c, c) for c in color_alt]
        
        cost_max = max(cost_acrossK.values())
        cost_min = min(cost_acrossK.values())
        
        for i_K, K in enumerate(K_range[:-1]):
            for i_mode in range(mode_numbers[i_K]):
                for i_mode2 in range(mode_numbers[i_K+1]):
                    m1, m2 = names_list[i_K][i_mode], names_list[i_K+1][i_mode2]
                    pair_label = "{}-{}".format(m1, m2)
                    if alt_color:
                        line_cmap = cmaps_alt[(i_mode+i_mode2+i_K+i_K+1) % len(cmaps_alt)]
                    elif line_cmap is None:
                        line_cmap = cm.get_cmap("Greys")

                    if pair_label in cost_acrossK:
                        ax = axes_handles[(K, i_mode)]
                        ax2 = axes_handles[(K_range[i_K+1], i_mode2)]
                        cost = float(cost_acrossK[pair_label])
                        edge_w = 0.85-(cost-cost_min)/(cost_max-cost_min)*0.8
                        
                        con = ConnectionPatch(xyA=(0.5,-0.05), coordsA='axes fraction', axesA=ax,
                                              xyB=(0.5,1.05), coordsB='axes fraction', axesB=ax2,
                                              color=line_cmap(edge_w), lw=4, alpha=1, zorder=0)
                        ax2.add_artist(con)
                        
                        # add text annotation  at midpoint
                        xyA_fig = fig.transFigure.inverted().transform(ax.transData.transform((N//2,-0.05)))
                        xyB_fig = fig.transFigure.inverted().transform(ax2.transData.transform((N//2,1.05)))
                        frac = 0.3 if (i_mode+i_mode2)%2==1 else 0.7 # to avoid overlapping texts
                        if mode_numbers[i_K]==1 or mode_numbers[i_K+1]==1:
                            frac = 0.5
                        mid_fig = get_frac_coord(xyA_fig,xyB_fig,frac=frac)
                        fig.text(mid_fig[0], mid_fig[1], "{:.3f}".format(cost), 
                                 fontsize=14, bbox=textbox_props, 
                                 ha='center', va='center', clip_on=False, zorder=1000)
                        
                        con.set_in_layout(False)

    return fig
            

def plot_alignment_list(mode_K: list[int], mode_names: list[str], cmap: list[tuple[float, float, float]], 
                   alignment_acrossK: dict, all_modes_alignment: dict, 
                   marker_size: float=200):
    """Plot alignment pattern (as a list) with all modes arranged in a column.
    Args:
        mode_K: List of K values for each mode.
        mode_names: List of mode names.
        cmap: Color map for clusters.
        alignment_acrossK: Alignment mappings across K.
        all_modes_alignment: Alignments within each mode.
        marker_size: Size of the markers in the plot. Defaults to 200.
    """
    K_max = np.max(mode_K)

    fig, ax = plt.subplots(1,1,figsize=(K_max*0.75,len(mode_names)*0.75), dpi=150)
    for i_m, mode_name in enumerate(mode_names[:-1]):
        K = mode_K[i_m]
        mode_pair = '{}-{}'.format(mode_name,mode_names[i_m+1])
        reordering_cur = all_modes_alignment[mode_name]
        ax.scatter(reordering_cur, np.ones(K)*i_m, s=marker_size, linewidths=0.5, edgecolors='k',
                   c=[cmap[reordering_cur[i]] for i in range(K)], zorder=4)
        mapping = alignment_acrossK[mode_pair]
        reordering_next = all_modes_alignment[mode_names[i_m+1]]
        for kp1 in range(len(mapping)):
            # ax.plot([reordering_cur[mapping[kp1]],reordering_next[kp1]], [i_m,i_m+1],
            #         c='k', ls='-', lw=0.8, zorder=5)
            ax.plot([reordering_cur.index(mapping[kp1]),reordering_next.index(kp1)], [i_m,i_m+1],
                    c='k', ls='-', lw=0.8, zorder=5)
    i_m += 1
    mode_name = mode_names[i_m]
    K = mode_K[i_m]
    reordering_cur = all_modes_alignment[mode_name]
    ax.scatter(reordering_cur, np.ones(K)*i_m, s=marker_size, linewidths=0.5, edgecolors='k',
            c=[cmap[reordering_cur[i]] for i in range(K)], zorder=4)

    ax.set_yticks(np.arange(len(mode_names)))
    ax.set_yticklabels(mode_names, rotation=0, fontsize=12, va='center', ha='right')
    ax.set_ylim(-0.5,len(mode_names)-0.5)
    ax.set_xlim(-0.5,K-0.5)
    ax.set_xticks(np.arange(K_max))
    ax.set_xticklabels([])
    ax.tick_params(axis='x',length=0)
    ax.set_xlabel("Clusters", fontsize=12)
    ax.set_ylabel("Modes", fontsize=12)
    ax.invert_yaxis()    

    return fig



def plot_alignment_graph(K_range: list[int], names_list: list[list[str]], cmap: list[tuple[float, float, float]], 
                   alignment_acrossK: dict, all_modes_alignment: dict, 
                   anchor_pairs: Optional[list[str]]=None,
                   wspace_padding: float=1.3, y_aspect: float=3, 
                   alt_color: bool=True,
                   color_alt: list[str]=['rosybrown', 'steelblue','goldenrod', 'darkseagreen'], 
                   ls_alt: list[str]=['-', '--'], 
                   marker_size: float=200, separate_labels: bool=False):
    """Plot alignment pattern (as a graph), with modes in each K arranged in a row, and space between modes.

    Args:
        K_range: List of K values.
        names_list: List of lists of mode names for each K.
        cmap: Color map for clusters.
        alignment_acrossK: Alignment mappings across K.
        all_modes_alignment: Global alignment (final output).
        wspace_padding: Horizontal space padding between modes. Defaults to 1.3.
        y_aspect: Aspect ratio for y-axis. Defaults to 3.
        color_alt: List of colors to alternate between for lines. Defaults to ['rosybrown', 'steelblue', 'goldenrod', 'darkseagreen'].
        ls_alt: List of line styles to distinguish lines between the anchor pair and the rest. Defaults to ['-', '--'].
        marker_size: float=200
    """

    mode_numbers = [len(names_list[i_K]) for i_K, K in enumerate(K_range)]
    n_row = len(K_range)
    n_col = max(mode_numbers)
    K_max = np.max(K_range)
    anchor_pairs_tuple = list()
    if anchor_pairs is not None:
        logger.info(f"Highlighting anchor pairs")
        for pr in anchor_pairs:
            M1, M2 = pr.split('-')
            anchor_pairs_tuple.append((M1, M2))

    # get the starting (x-)positions of each mode in each row
    start_positions = np.zeros((n_row, n_col), dtype=int)
    for i_K, K in enumerate(K_range):
        for i_mode in range(n_col):
            start_positions[i_K, i_mode] = i_mode * int(K_max*wspace_padding)  # add some space between modes

    fig, ax = plt.subplots(1,1,figsize=(n_col*K_max*0.3,n_row*1.5), dpi=150)
    for i_K, K in enumerate(K_range[:-1]):
        i_K2, K2 = i_K+1, K_range[i_K+1]
        for i_mode in range(mode_numbers[i_K]):
            mode_name = names_list[i_K][i_mode]
            ax.scatter(all_modes_alignment[mode_name]+start_positions[i_K, i_mode], np.ones(K)*i_K, 
                       s=marker_size, linewidths=0.5, edgecolors='k',
                    c=[cmap[all_modes_alignment[mode_name][i]] for i in range(K)], zorder=6)
            for i_mode2 in range(mode_numbers[i_K2]):
                mode_name2 = names_list[i_K2][i_mode2]
                mode_pair = '{}-{}'.format(mode_name,mode_name2)
                if mode_pair in alignment_acrossK:
                    mapping = alignment_acrossK[mode_pair]
                    reordering_cur = all_modes_alignment[mode_name]
                    reordering_next = all_modes_alignment[mode_name2]
                    for kp1 in range(len(mapping)):
                        # switch between colors for better visibility
                        if alt_color:
                            color = color_alt[(i_mode+i_mode2+i_K+i_K2) % len(color_alt)]
                        else:
                            color = 'k'
                        # cycle between line styles for better visibility
                        if anchor_pairs is not None and (mode_name, mode_name2) in anchor_pairs_tuple:
                            ls = ls_alt[0]
                        else:
                            ls = ls_alt[1] #ls_alt[(i_mode+i_mode2) % len(ls_alt)]
                        if reordering_cur.index(mapping[kp1])!=reordering_next.index(kp1):
                            ax.plot([reordering_cur.index(mapping[kp1])+start_positions[i_K, i_mode],
                                     reordering_next.index(kp1)+start_positions[i_K2, i_mode2]], 
                                     [i_K+0.1,i_K2-0.1],
                                    c=color, ls=ls, lw=0.8, zorder=2)
                        else:
                            ax.plot([reordering_cur.index(mapping[kp1])+start_positions[i_K, i_mode],
                                     reordering_next.index(kp1)+start_positions[i_K2, i_mode2]], 
                                     [i_K+0.1,i_K2-0.1],
                                    c='lightgrey', ls=':', lw=0.3, zorder=2)
                    # add a rectangle to highlight the mode
                    if separate_labels:
                        ax.text(start_positions[i_K2, i_mode2], i_K2+0.2, mode_name2.replace('_',' '),
                                c='gray', fontsize=9, ha='left', va='center')
                    rect = Rectangle((start_positions[i_K2, i_mode2]-0.5, i_K2-0.15), K2, 0.3, 
                                    linewidth=0.5, edgecolor='lightgrey', facecolor='lightgrey', alpha=0.1, 
                                    joinstyle='round', capstyle='round',zorder=1)
                    ax.add_patch(rect)
                else:
                    print(f"Warning: alignment for mode pair {mode_pair} not found.")

    i_K = len(K_range)-1
    K = K_range[i_K]
    for i_mode in range(mode_numbers[i_K]):
        mode_name = names_list[i_K][i_mode]
        ax.scatter(all_modes_alignment[mode_name]+start_positions[i_K, i_mode], np.ones(K)*i_K, 
                   s=marker_size, linewidths=0.5, edgecolors='k',
                    c=[cmap[all_modes_alignment[mode_name][i]] for i in range(K)], zorder=6)
    i_K = 0
    K = K_range[i_K]
    for i_mode in range(mode_numbers[i_K]):
        mode_name = names_list[i_K][i_mode]   
        if separate_labels:
            ax.text(start_positions[i_K, i_mode], i_K+0.2, mode_name.replace('_',' '),
                    c='gray', fontsize=9, ha='left', va='center')   
        # add a rectangle to highlight the mode
        rect = Rectangle((start_positions[i_K, i_mode]-0.5, i_K-0.15), K, 0.3, 
                        linewidth=0.5, edgecolor='lightgrey', facecolor='lightgrey', alpha=0.1, 
                        joinstyle='round', capstyle='round',zorder=1)
        ax.add_patch(rect)

    tick_positions = [start_positions[0, i_mode] + (K_max - 1) / 2.0 for i_mode in range(n_col)]
    ax.set_xticks(tick_positions)
    if not separate_labels:
        ax.set_xticklabels([f"M{i_m+1}" for i_m in range(n_col)], fontsize=14)
    ax.set_xlabel("Clusters in Modes", fontsize=14)
    ax.set_xlim(-1, (n_col-1) * int(K_max*wspace_padding) + K_max )  

    ax.set_yticks(np.arange(n_row))
    if not separate_labels:
        ax.set_yticklabels([f"K{K}" for K in K_range], rotation=0, fontsize=14, va='center', ha='right')
    ax.set_ylim(-0.5,n_row-0.5)

    ax.invert_yaxis() 
    ax.set_aspect(y_aspect, adjustable='box')
    return fig