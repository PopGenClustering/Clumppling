
import numpy as np
import os
import matplotlib.pyplot as plt
from itertools import product,combinations_with_replacement
from collections import defaultdict
import networkx as nx

def alignQ_wrtP(P,Q,idxQ2P,merge=True):
    # K1, K2 = P.shape[1], Q.shape[1]
    idxQ2P = list(idxQ2P)
    if merge:        
        aligned_Q = np.zeros_like(P)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        dups = np.unique([i for i in idxQ2P if idxQ2P.count(i)>1])
        extras = list()
        dups_min = defaultdict(lambda: (float('inf'),None))

        new_pattern = [0 for _ in range(Q.shape[1])]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx not in dups:
                new_pattern[q_idx] = p_idx
                aligned_Q[:,p_idx] = Q[:,q_idx]
            else:
                diff = np.linalg.norm(Q[:,q_idx]-P[:,p_idx])
                if dups_min[p_idx][0] > diff:
                    dups_min[p_idx] = (diff,q_idx) 
        extra_cnt = P.shape[1]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx in dups:
                if q_idx==dups_min[p_idx][1]:
                    new_pattern[q_idx] = p_idx
                    aligned_Q[:,p_idx] = Q[:,q_idx]
                else:
                    new_pattern[q_idx] = extra_cnt
                    extras.append(q_idx)
                    extra_cnt += 1
        
        for ie, e in enumerate(extras):
            # print(P.shape[1],ie,e)
            # print(aligned_Q.shape)
            # print(P.shape)
            # print(Q.shape)
            aligned_Q[:,P.shape[1]+ie] = Q[:,e]
            
    return aligned_Q, new_pattern


def plot_membership(ax,P,max_K,cmap,title):

    N = P.shape[0]
    P_aug = np.hstack((np.zeros((N,1)),P))

    K = P.shape[1]
    for k in range(K):
        ax.bar(range(N), P_aug[:,(k+1)], bottom=np.sum(P_aug[:,:(k+1)],axis=1), 
               width=1.0, edgecolor='w', linewidth=0.05, facecolor=cmap(k/max_K))

    ax.set_xticks([])
    ax.set_xlim([0,N])
    ax.set_ylim([0,1])
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([])
    # ax.set_title(title)
    if title:
        ax.set_ylabel("\n".join(title.split()), rotation=0, fontsize=18, labelpad=30, va="center" )
    else:
        ax.set_ylabel("")
    return

def plot_aligned(K,m,Q_list,modes,align_ILP_res,rep_modes,consensusQ,max_K,k2ids,idx2idxinK,save_path,plot_name,cmap):
    # consensusQ = consensusQ_modes[K][m]
    consR_idx = rep_modes[K][m]
    mode_size = len(modes[m])
    fig, axes = plt.subplots(mode_size+1,1,figsize=(20,2*(mode_size+1)))
    ax = axes[0]
   
    plot_membership(ax,consensusQ,max_K,cmap,"K{} mode{} average membership".format(K,m))
    for i,r in enumerate(modes[m]):
        ax = axes[i+1]
        Q = Q_list[k2ids[K][0]+r]
        if r!=idx2idxinK[consR_idx]:
            aligned_idxQ2P = align_ILP_res[K][r][idx2idxinK[consR_idx]]
            aligned_Q = np.zeros_like(consensusQ)
            for q_idx in range(Q.shape[1]):
                aligned_Q[:,aligned_idxQ2P[q_idx]] += Q[:,q_idx]
        else:
            aligned_Q = Q
        if i==0:
            title = "K{} mode{} replicates".format(K,m)
        else:
            title = None
        plot_membership(ax,aligned_Q,max_K,cmap,title)        
    # plt.suptitle("K{} mode{} Average Membership".format(K,m))
    fig.savefig(os.path.join(save_path,plot_name),dpi=50)
    plt.close(fig)


def plot_acrossK_multipartite(K_range,modes_allK_list,meanQ_modes,meanQ_acrossK_cost,layer_color,title,save_path,plot_file_name):
    
    K_range_sorted = sorted(K_range,reverse=True)
    max_cost = max([float(c) for c in meanQ_acrossK_cost.values()])
    max_mode_num = 0
    mode_sizes = defaultdict(list)
    
    G = nx.Graph()
    # add nodes 
    nid = 0
    km2nid = dict()
    
    for K in K_range_sorted:
        for m in range(len(modes_allK_list[K])):
            G.add_node(nid, layer=K, label="K{}#{}".format(K,m))
            km2nid["{}#{}".format(K,m)] = nid
            nid += 1
            mode_sizes[K].append(len(modes_allK_list[K][m]))
        if len(meanQ_modes[K])>max_mode_num:
            max_mode_num = len(meanQ_modes[K])
    
    nid = 0
    nodes_sizes = list()
    for K in K_range_sorted:
        for m in range(len(modes_allK_list[K])):   
            s = mode_sizes[K][m]/sum(mode_sizes[K])
            G.nodes[nid]["size"] = s
            nodes_sizes.append(s)
            nid += 1
    
    # add edges (edge weight: 1-norm cost)
    for k1k2_pair in meanQ_acrossK_cost:
        k1m = k1k2_pair.split("-")[0]
        k2m = k1k2_pair.split("-")[1]
        if int(k1m.split("#")[0])+1==int(k2m.split("#")[0]):
            cost = float(meanQ_acrossK_cost[k1k2_pair])
            norm_cost = 1-cost/max_cost
            G.add_edge(km2nid[k1m],km2nid[k2m],weight=norm_cost,cost=round(cost,3))
    
    # plot the network of modes 
    color = [layer_color[data["layer"]] for v, data in G.nodes(data=True)]
    size = [s * 2500 for s in nodes_sizes]
    fig, ax = plt.subplots(figsize=(2*len(K_range_sorted), 2*max_mode_num))
    
    pos = nx.multipartite_layout(G, subset_key="layer")
    node_labels = nx.get_node_attributes(G, 'label') 
    weights = list(nx.get_edge_attributes(G,'weight').values())
    edge_labels = nx.get_edge_attributes(G,'cost')
    nx.draw(G, pos, node_color=color, node_size=size, 
            font_size=12, labels=node_labels, with_labels = True,
            width=2, edge_color=weights, edge_cmap=plt.cm.Reds)
    
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, label_pos=0.73, alpha=0.7, font_size=10)
    # plt.title(title, fontsize=16)
    # save 
    fig.savefig(os.path.join(save_path,plot_file_name))
    plt.close(fig)
    
    return G


def plot_acrossK_chains(K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_acrossK_cost,save_path,plot_file_name_suffix,cmap,merge_cluster=False):
    K_range_sorted = sorted(K_range,reverse=True)
    best_alignments_cost = dict()
    best_alignments = dict()
    for k1k2_pair in meanQ_acrossK_cost:
        k1m = k1k2_pair.split("-")[0]
        k2m = k1k2_pair.split("-")[1]
        k1 = int(k1m.split("#")[0])
        k2 = int(k2m.split("#")[0]) # k1<k2
        if K_range_sorted.index(k1)-1==K_range_sorted.index(k2):
            cost = float(meanQ_acrossK_cost[k1k2_pair])
            if k1m not in best_alignments_cost:
                best_alignments_cost[k1m] = cost
                best_alignments[k1m] = k2m
            else:
                if cost<best_alignments_cost[k1m]:
                    best_alignments_cost[k1m] = cost
                    best_alignments[k1m] = k2m
    
    # retrieve optimally aligned modes in chain for each mode of min(K)
    min_K, max_K = min(K_range), max(K_range)
    best_alignment_chains = dict()
    for minK_m in meanQ_modes[min_K]:
        mode_lb = "{}#{}".format(min_K,minK_m)
        starting_mode = mode_lb
        chain = [mode_lb]
        while mode_lb in best_alignments:
            mode_lb = best_alignments[mode_lb]
            chain.append(mode_lb)
        best_alignment_chains[starting_mode] = chain
    
    # retrieve optimal alignment
    for minK_m in range(len(meanQ_modes[min_K].keys())):
        mode_lb = "{}#{}".format(min_K,minK_m)
        chain = best_alignment_chains[mode_lb]
        alignment_pattern = list()
        for i_c in range(len(chain)-1):
            m1,m2 = chain[i_c], chain[i_c+1]
            if i_c==0:
                pattern = meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]
            else:
                prev_pattern = alignment_pattern[i_c-1]
                pattern = [prev_pattern[i] for i in meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]]
            if not merge_cluster:
                pattern_expanded = list()
                extra_cnt = len(np.unique(pattern))
                for p in pattern:
                    if p not in pattern_expanded:
                        pattern_expanded.append(p)
                    else:
                        pattern_expanded.append(extra_cnt)
                        extra_cnt += 1
                alignment_pattern.append(pattern_expanded)
            else:
                alignment_pattern.append(pattern)
        
        # plot chains of modes
        fig, axes = plt.subplots(len(chain),1,figsize=(20,2*len(chain)))
        ax = axes[0]
        K,m = int(chain[0].split("#")[0]),int(chain[0].split("#")[1])
        plot_name = "K{}_mode{}_{}.png".format(K,m,plot_file_name_suffix)
        Q = meanQ_modes[K][m]
        plot_membership(ax,Q,max_K,cmap,"K{} mode{}".format(K,m))
        
        
        for i,ali_pat in enumerate(alignment_pattern):
            ax = axes[i+1]
            m_next = chain[i+1]
            K,m = int(m_next.split("#")[0]),int(m_next.split("#")[1])
            Q = meanQ_modes[K][m]
            aligned_Q = np.zeros_like(Q)
            for q_idx in range(Q.shape[1]):
                aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
            plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
        
        fig.savefig(os.path.join(save_path,plot_name),dpi=50)
        plt.close(fig)
    
    return best_alignment_chains    
    

def plot_withinK_modes(K,max_K,meanQ_modes,meanQ_acrossK_Q2P,save_path,plot_file_name_suffix,cmap):
    
    modes = range(len(meanQ_modes[K].keys()))
    fig, axes = plt.subplots(len(modes),1,figsize=(20,2*len(modes)))
    if len(modes)==1:
        ax = axes
    else:
        ax = axes[0]
    plot_name = "K{}_{}.png".format(K,plot_file_name_suffix)
    m = 0
    m1 = "{}#{}".format(K,m)
    Q = meanQ_modes[K][m]
    plot_membership(ax,Q,max_K,cmap,"K{} mode{}".format(K,m))
    
    
    for i in range(1,len(modes)):
        m = modes[i]
        m2 = "{}#{}".format(K,m)
        # retrieve alignment pattern
        pattern = meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]
        Q = meanQ_modes[K][m]
        aligned_Q = np.zeros_like(Q)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pattern[q_idx]] += Q[:,q_idx]
        # plot
        ax = axes[i]
        plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
    
    fig.savefig(os.path.join(save_path,plot_name),dpi=50)
    plt.close(fig)
    
    
def show_all_modes_old(plot_flag_all_modes,K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_best_ILP_acrossK,save_path,plot_name,cmap):
    
    # use mode 0 as base for each K
    # align all mode 0 across K's
    # use the best alignment mode as the ref
    meanQ_best_ILP_acrossK = meanQ_best_ILP_acrossK[::-1]
    num_modes = sum([len(meanQ_modes[K].keys()) for K in K_range])  
    mode2fig_idx = ["{}#{}".format(K,m) for K in K_range for m in meanQ_modes[K].keys()]
    mode2fig_idx = {m:i for i,m in enumerate(mode2fig_idx)}
    
    modeQ_path = os.path.join(save_path,"modes_Q")
    
    if plot_flag_all_modes:
        fig, axes = plt.subplots(num_modes,1,figsize=(15,1.5*num_modes))
    
    K = K_range[0]
    # m = 0
    m1 = meanQ_best_ILP_acrossK[0].split("-")[0]
    assert(m1.split("#")[0]==str(K))
    m_m1 = int(m1.split("#")[1])
    # m1 = "{}#{}".format(K,m)
    max_K = max(K_range)
    base_Q = meanQ_modes[K][m_m1]
    
    if plot_flag_all_modes:
        ax = axes[mode2fig_idx[m1]]
        plot_membership(ax,base_Q,max_K,cmap,"K{} mode{}".format(K,m_m1)) 
    np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m_m1)), base_Q, delimiter=' ')
    
    # fig_idx = 1
    modes = range(len(meanQ_modes[K].keys()))
    if len(modes)>1:
        for m in modes:
            if m!=m_m1:
                m2 = "{}#{}".format(K,m)
                # retrieve alignment pattern
                ali_pat = meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]
                Q = meanQ_modes[K][m]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                if plot_flag_all_modes:
                    ax = axes[mode2fig_idx[m2]]#axes[fig_idx]
                    plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
                np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m)), aligned_Q, delimiter=' ')
                # fig_idx += 1
    
    prev_pattern = None
    prev_mode = None
    # K_idx = 1
    # while K_idx<len(K_range):
    for l in meanQ_best_ILP_acrossK:
        
        m1 = l.split("-")[0]#"{}#{}".format(K_range[K_idx-1],0)
        m2 = l.split("-")[1]
        
        K = int(m2.split("#")[0])
        
        modes = range(len(meanQ_modes[K].keys()))
        m_m2 = int(m2.split("#")[1])
        Q = meanQ_modes[K][m_m2]
        # m2 = "{}#{}".format(K,m)
        if prev_mode and m1!=prev_mode:
            pattern = meanQ_acrossK_Q2P["{}-{}".format(prev_mode,m1)]
            if prev_pattern:
                prev_pattern = [prev_pattern[i] for i in pattern]
            else:
                prev_pattern = pattern
        # retrieve alignment pattern
        pattern = meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]
        
        
        # P = meanQ_modes[int(m1.split("#")[0])][int(m1.split("#")[1])]
        # # if prev_pattern:
        # #     P = P[:,prev_pattern]
        #     # reord = [prev_pattern[i] for i in pattern]
        #     # pattern = [prev_pattern[i] for i in pattern]
        # aligned_Q, new_pattern = alignQ_wrtP(P,Q,pattern,merge=False)
        # if prev_pattern:
        #     new_pattern = [prev_pattern[i] for i in new_pattern]
        
 
        if prev_pattern:
            pattern = [prev_pattern[i] for i in pattern]
            
        
        pattern_expanded = list()
        extra_cnt = len(np.unique(pattern))
        for p in pattern:
            if p not in pattern_expanded:
                pattern_expanded.append(p)
            else:
                pattern_expanded.append(extra_cnt)
                extra_cnt += 1
        
        aligned_Q = np.zeros_like(Q)
        
        
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pattern_expanded[q_idx]] += Q[:,q_idx]
            
        if plot_flag_all_modes:
            ax = axes[mode2fig_idx[m2]] #axes[fig_idx]        
            plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m_m2))
        np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m_m2)), aligned_Q, delimiter=' ')
        # fig_idx += 1
        prev_pattern = pattern_expanded  
        
        
        if len(modes)>1:
            for m in modes:
                if m!=m_m2:
                    m3 = "{}#{}".format(K,m)
                    ali_pat = meanQ_acrossK_Q2P["{}-{}".format(m2,m3)]
                    ali_pat = [prev_pattern[i] for i in ali_pat]
                    Q = meanQ_modes[K][m]
                    aligned_Q = np.zeros_like(Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                    if plot_flag_all_modes:
                        ax = axes[mode2fig_idx[m3]]
                        # ax = axes[fig_idx]    
                        plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
                        # fig_idx += 1
                    np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m)), aligned_Q, delimiter=' ')            
                        
        prev_mode = m2      
                     
    if plot_flag_all_modes:    
        fig.savefig(os.path.join(save_path,plot_name),dpi=50)
        plt.close(fig)
        
    # return
    
    
def show_all_modes(plot_flag_all_modes,K_range,meanQ_modes,meanQ_acrossK_Q2P,meanQ_best_ILP_acrossK,save_path,plot_name,cmap):
    
    # use mode 0 as base for each K
    # align all mode 0 across K's
    # use the best alignment mode as the ref
    meanQ_best_ILP_acrossK = meanQ_best_ILP_acrossK[::-1]
    num_modes = sum([len(meanQ_modes[K].keys()) for K in K_range])  
    mode2fig_idx = ["{}#{}".format(K,m) for K in K_range for m in meanQ_modes[K].keys()]
    mode2fig_idx = {m:i for i,m in enumerate(mode2fig_idx)}
    
    modeQ_path = os.path.join(save_path,"modes_Q")
    base_patterns = dict()
    
    if plot_flag_all_modes:
        fig, axes = plt.subplots(num_modes,1,figsize=(15,1.5*num_modes))
    
    K = K_range[0]
    # m = 0
    m1 = meanQ_best_ILP_acrossK[0].split("-")[0]
    assert(m1.split("#")[0]==str(K))
    m_m1 = int(m1.split("#")[1])
    # m1 = "{}#{}".format(K,m)
    max_K = max(K_range)
    base_Q = meanQ_modes[K][m_m1]
    base_patterns[m1] = [i for i in range(K)]
    if plot_flag_all_modes:
        ax = axes[mode2fig_idx[m1]]
        plot_membership(ax,base_Q,max_K,cmap,"K{} mode{}".format(K,m_m1)) 
    np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m_m1)), base_Q, delimiter=' ')
    
    modes = range(len(meanQ_modes[K].keys()))
    if len(modes)>1:
        for m in modes:
            if m!=m_m1:
                m2 = "{}#{}".format(K,m)
                # retrieve alignment pattern
                ali_pat = meanQ_acrossK_Q2P["{}-{}".format(m1,m2)]
                Q = meanQ_modes[K][m]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                base_patterns[m2] = ali_pat
                if plot_flag_all_modes:
                    ax = axes[mode2fig_idx[m2]]#axes[fig_idx]
                    plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
                np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m)), aligned_Q, delimiter=' ')
    
    for best_pair in meanQ_best_ILP_acrossK:
        m1 = best_pair.split("-")[0]#"{}#{}".format(K_range[K_idx-1],0)
        m2 = best_pair.split("-")[1]
        
        K = int(m2.split("#")[0])
        
        modes = range(len(meanQ_modes[K].keys()))
        m_m2 = int(m2.split("#")[1])
        Q = meanQ_modes[K][m_m2]
        P = meanQ_modes[int(m1.split("#")[0])][int(m1.split("#")[1])]
        pattern = meanQ_acrossK_Q2P[best_pair]
        aligned_Q, new_pattern = alignQ_wrtP(P,Q,pattern,merge=False)

        pat = [base_patterns[m1][i] if i<int(m1.split("#")[0]) else i for i in new_pattern]
        aligned_Q = np.zeros_like(aligned_Q)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
        base_patterns[m2] = pat
        
        if plot_flag_all_modes:
            ax = axes[mode2fig_idx[m2]] #axes[fig_idx]        
            plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m_m2))
        np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m_m2)), aligned_Q, delimiter=' ')
        
        if len(modes)>1:
            for m in modes:
                if m!=m_m2:
                    m3 = "{}#{}".format(K,m)
                    ali_pat = meanQ_acrossK_Q2P["{}-{}".format(m2,m3)]
                    Q = meanQ_modes[K][m]
                    P = meanQ_modes[K][m_m2]
                    aligned_Q, new_pattern = alignQ_wrtP(P,Q,ali_pat,merge=False)
                    
                    pat = [base_patterns[m2][i] for i in new_pattern ]
                    aligned_Q = np.zeros_like(aligned_Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
                    base_patterns[m3] = pat
                    
                    if plot_flag_all_modes:
                        ax = axes[mode2fig_idx[m3]]
                        # ax = axes[fig_idx]    
                        plot_membership(ax,aligned_Q,max_K,cmap,"K{} mode{}".format(K,m))
                        # fig_idx += 1
                    np.savetxt(os.path.join(modeQ_path,'K{}_mode{}.Q'.format(K,m)), aligned_Q, delimiter=' ')            
        
    if plot_flag_all_modes:    
        fig.savefig(os.path.join(save_path,plot_name),dpi=50)
        plt.close(fig)
    

                     
    