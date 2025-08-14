import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
# import matplotlib.colors as mcolors
from matplotlib.patches import ConnectionPatch
from typing import Optional, cast
from clumppling.utils import get_uniq_lb_sep

import logging
logger = logging.getLogger(__name__)

def plot_colorbar(rgb_colors: list[tuple[float, float, float]], K_max: int, fig_dir:str):
    assert K_max == len(rgb_colors), f"Number of colors in colormap ({len(rgb_colors)}) is K_max ({K_max})"

    fig, ax = plt.subplots(figsize=(K_max, 1),facecolor='white',subplot_kw=dict(yticks=[])) 
    ax.set_xticks(ticks=np.arange(0.5,K_max+0.5,1),labels=["Cls.{}".format(k) for k in np.arange(1,K_max+1)])
    img = np.array([rgb_colors])
    ax.imshow(img, aspect='auto', extent=(0, K_max, 0, 1))
    ax.get_yaxis().set_visible(False)

    fig.savefig(os.path.join(fig_dir,"colorbar.png"), bbox_inches='tight',dpi=300, transparent=False)
    plt.close(fig)


# plot_modes_withinK,plot_modes, plot_major_modes, plot_all_modes,

def plot_membership(Q:np.ndarray, cmap: list[tuple[float, float, float]], 
                    ax: Optional[maxes.Axes] = None, ylab: str="", title: str="",
                    fontsize: float=14) -> maxes.Axes:
    if ax is None:
        _, ax = plt.subplots()
        ax = cast(maxes.Axes, ax)

    N = Q.shape[0]
    Q_aug = np.hstack((np.zeros((N,1)),Q))

    K = Q.shape[1]
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


def plot_membership_reordered(Q:np.ndarray, lbs: list, cmap: list[tuple[float, float, float]], 
                              ax: Optional[maxes.Axes] = None, ylab: str="", title: str="", 
                              fontsize: float=14) -> maxes.Axes:
    if ax is None:
        _, ax = plt.subplots()
        ax = cast(maxes.Axes, ax)
    N = Q.shape[0]
    K = Q.shape[1]
    mbsp_sum = Q.sum(axis=0)
    mbsp_sortidx = np.argsort(-mbsp_sum) # from largest to smallest
    Q = Q[:,mbsp_sortidx]

    for lb in np.unique(lbs):
        lb_indices = np.where(lbs==lb)[0]
        lb_P = Q[lb_indices,:]
        largest_cls_mbsp = lb_P[:,0]
        ind_sortidx = np.argsort(-largest_cls_mbsp)
        lb_P_sorted = lb_P[ind_sortidx,:]
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
                          names: list[str], ind_labels: list=[]):
    # get separation points in the labels
    if ind_labels:
        uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx = get_uniq_lb_sep(ind_labels)

    fig, axes = plt.subplots(len(Q_list),1,figsize=(15,1.5*len(Q_list)),facecolor='white')

    for i in range(len(Q_list)):
        Q = Q_list[i]
        lb = names[i]
        ax = axes[i] if len(Q_list)>1 else axes
        plot_membership(Q, cmap, ax=ax, ylab=lb)
        if ind_labels:
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
    
    return fig


def plot_graph(K_range: list[int], Q_list_list: list[list[np.ndarray]], cmap: list[tuple[float, float, float]], 
               names_list: list[list[str]], labels_list: Optional[list[list[str]]]=None,
               cost_acrossK: Optional[dict]=None,
               ind_labels: list=[], fontsize: float=14, line_cmap=None):
    mode_numbers = [len(names_list[i_K]) for i_K, K in enumerate(K_range)]
    n_row = len(K_range)
    n_col = max(mode_numbers)
    N = Q_list_list[0][0].shape[0]

    # get separation points in the labels
    if ind_labels:
        uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx = get_uniq_lb_sep(ind_labels)
    if labels_list is None:
        labels_list = names_list

    fig = plt.figure(figsize=(7*n_col,2*n_row), dpi=150)
    gs = fig.add_gridspec(n_row,n_col, wspace=0.2, hspace=1)

    axes_handles = dict()
    # modes_handles = dict()
    for i_K, K in enumerate(K_range):
        for i_mode in range(mode_numbers[i_K]):
            ax = fig.add_subplot(gs[i_K, i_mode])
            axes_handles[(K, i_mode)] = ax
            # modes_handles[names_list[i_K][i_mode]] = ax
            Q_label = labels_list[i_K][i_mode] 
            Q = Q_list_list[i_K][i_mode]
            plot_membership(Q, cmap, ax=ax, title=Q_label, fontsize=fontsize)
            ax.set_zorder(9)
            if ind_labels:
                if len(ind_labels)!=N:
                    logger.warning(f"Number of individual labels ({len(ind_labels)}) does not match number of rows in Q matrix ({N}). ")
                for v in uniq_lbs_sep_idx:
                    ax.axvline(v, ymin=-0.2, ymax=1, color='k', ls='--', lw=0.5, clip_on=False)
                if i_K==len(K_range)-1:
                    ax.set_xticks(uniq_lbs_indices)
                    ax.tick_params(axis='x', which='both', length=0)
                    rot = 45 if np.any([len(lb)>5 for lb in uniq_lbs]) else 0
                    fs = 14 if len(uniq_lbs)<10 else 10
                    ax.set_xticklabels(uniq_lbs, rotation=rot, ha='center', fontsize=fs)

    for i_K, K in enumerate(K_range):
        ax = axes_handles[(K, 0)]
        ax.set_ylabel(f"K={K}", fontsize=18, rotation=0, labelpad=10, va="center", ha="right", weight='bold')

    # plot cost if applicable
    def get_frac_coord(coord1,coord2,frac: float=0.5):
        return (coord2[0]-coord1[0])*frac+coord1[0], (coord2[1]-coord1[1])*frac+coord1[1] 
    
    if cost_acrossK is not None:
        
        logger.info(f"Plot cost across K as graph edges")
        textbox_props = dict(boxstyle='round', facecolor='white', alpha=0.6, edgecolor='none', pad=0.1)
        
        if line_cmap is None:
            logger.info(f"Using default Greys colormap for cost lines")
            line_cmap = plt.get_cmap("Greys")
        cost_max = max(cost_acrossK.values())
        cost_min = min(cost_acrossK.values())
        
        for i_K, K in enumerate(K_range[:-1]):
            for i_mode in range(mode_numbers[i_K]):
                for i_mode2 in range(mode_numbers[i_K+1]):
                    m1, m2 = names_list[i_K][i_mode], names_list[i_K+1][i_mode2]
                    pair_label = "{}-{}".format(m1, m2)
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
                        frac = 0.3 if (i_mode)%2==1 else 0.7 # to avoid overlapping texts
                        if mode_numbers[i_K]==1 or mode_numbers[i_K+1]==1:
                            frac = 0.5
                        mid_fig = get_frac_coord(xyA_fig,xyB_fig,frac=frac)
                        fig.text(mid_fig[0], mid_fig[1], "{:.3f}".format(cost), 
                                 fontsize=14, bbox=textbox_props, 
                                 ha='center', va='center', clip_on=False, zorder=1000)
                        
                        con.set_in_layout(False)

    return fig
            

def plot_alignment(mode_K: list[int], mode_names: list[str], cmap: list[tuple[float, float, float]], 
                   alignment_acrossK: dict, all_modes_alignment: dict, 
                   marker_size: float=150):
    K_max = np.max(mode_K)

    if len(mode_names)>4:
        fig, ax = plt.subplots(1,1,figsize=(K_max//2,len(mode_names)), dpi=150)
    else:
        fig, ax = plt.subplots(1,1,figsize=(K_max//2,len(mode_names)*0.75), dpi=150)
    for i_m, mode_name in enumerate(mode_names[:-1]):
        K = mode_K[i_m]
        mode_pair = '{}-{}'.format(mode_name,mode_names[i_m+1])
        reordering_cur = all_modes_alignment[mode_name]
        ax.scatter(reordering_cur, np.ones(K)*i_m, s=marker_size, linewidths=0.5, edgecolors='k',
                   c=[cmap[reordering_cur[i]] for i in range(K)], zorder=4)
        mapping = alignment_acrossK[mode_pair]
        reordering_next = all_modes_alignment[mode_names[i_m+1]]
        for kp1 in range(len(mapping)):
            ax.plot([reordering_cur[mapping[kp1]],reordering_next[kp1]], [i_m,i_m+1],
                    c='k', ls='-', lw=0.8, zorder=5)
    i_m += 1
    mode_name = mode_names[i_m]
    K = mode_K[i_m]
    reordering_cur = all_modes_alignment[mode_name]
    ax.scatter(reordering_cur, np.ones(K)*i_m, s=marker_size, linewidths=0.5, edgecolors='k',
            c=[cmap[reordering_cur[i]] for i in range(K)], zorder=4)

    ax.set_yticks(np.arange(len(mode_names)))
    ax.set_yticklabels(mode_names, rotation=0, fontsize=12, va='center', ha='right')
    ax.set_ylim([-0.5,len(mode_names)-0.5])
    ax.set_xlim([-0.5,K-0.5])
    ax.set_xticks(np.arange(K_max))
    ax.set_xticklabels([])
    ax.tick_params(axis='x',length=0)
    ax.set_xlabel("Clusters", fontsize=12)
    ax.set_ylabel("Modes", fontsize=12)
    ax.invert_yaxis()    

    return fig