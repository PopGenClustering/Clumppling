import os
import numpy as np
import pandas as pd
import shutil
from scipy.spatial.distance import cdist
import cvxpy as cp
from collections import defaultdict, Counter
import networkx as nx
from cdlib import algorithms
from itertools import combinations, product 

from typing import List, Tuple, Any, TypedDict

from clumppling.utils import cost_membership, cost_membership_sep, pattern_to_str, disp_msg, C2Gprime, avg_tril
from .parseInput import validate_membership_matrix
from .alignWithinK import load_withinK_qfiles, write_aligned_within_k, aligned_res_to_dicts
from .detectMode import construct_cost_mat, community_labels_to_modes, cost_to_adj, test_comm_struc, ModesDict, cd_custom, compute_mode_avg_stats #, detect_communities, community_labels_to_modes

import logging
logger = logging.getLogger(__name__)


def align_ILP(P: np.ndarray, Q: np.ndarray) -> Tuple[Any, List[int]]:
    """ Align two matrices P and Q using Integer Linear Programming (ILP).
        P and Q are matrices with the same number of rows (individuals) but different number of columns (clusters).
        The goal is to find a mapping from clusters in Q to clusters in P that minimizes the squared Euclidean distance between the clusters.
        Args:
            P: Membership matrix of shape (N_ind, K_P) for clusters in P.
            Q: Membership matrix of shape (N_ind, K_Q) for clusters in Q.
        Returns:
            opt_obj: The optimal objective value (minimum distance).
            idxQ2P: List of indices mapping clusters in Q to clusters in P.
                    (e.g., idxQ2P = [1, 2, 0] means Q[0] -> P[1], Q[1] -> P[2], Q[2] -> P[0])
    """
    # Check dimensions
    assert P.ndim == 2 and Q.ndim == 2, "Both P and Q must be 2D matrices."
    assert P.shape[0] == Q.shape[0], "P and Q must have the same number of rows (individuals)."
    assert P.shape[1] <= Q.shape[1], "Q must have no fewer rows (clusters) than P."
    
    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    N_ind = Q.shape[0] #individuals
    
    # Compute distance
    D = cdist(Q.T, P.T, metric='sqeuclidean')/N_ind # seuclidean, euclidean  
    d = D.flatten()
    
    # Set constraints
    # sum of row = 1
    A = np.zeros((K_Q,d.shape[0]))
    for i in range(K_Q):
        A[i,i*K_P:(i+1)*K_P] = np.ones((1,K_P))
    A_ones = np.ones((K_Q))
    # sum of column >= 1
    B = np.zeros((K_P,d.shape[0]))
    for i in range(K_P):
        B[i,np.arange(i,d.shape[0],K_P)] = np.ones((1,K_Q))
    B_ones = np.ones((K_P))
 
    # Solve ILP
    w = cp.Variable(d.shape[0], boolean = True, name='w')
    constraints = [A @ w == A_ones, B @ w >= B_ones]
    obj = d.T @ w
    problem = cp.Problem(cp.Minimize(obj), constraints)
    problem.solve(solver=cp.GLPK_MI) #
    
    # Get optimal solution
    if problem.value is None:
        logger.error("ERROR solving the ILP. Check the input matrices P and Q.")
        return None, []
    else:
        opt_obj = problem.value

    var_dict = {var.name(): var for var in problem.variables()}
    opt_w = np.asarray(var_dict['w'].value, dtype=int)
    opt_W = np.reshape(opt_w, D.shape)

    # Match index Q_idx to P_idx
    idxQ2P = np.where(opt_W==1)[1]
    
    return opt_obj, list(idxQ2P)


# align P and Q by enumerating all combinations of two clusters when K is differ by 1
def align_ILP_diff1(P: np.ndarray, Q: np.ndarray) -> Tuple[Any, List[int]]:

    K_Q = Q.shape[1] # #rows
    K_P = P.shape[1] # #columns
    assert K_Q-K_P==1, "K of Q must be exactly one more than K of P."

    # all combinations of Q into the same size of P
    Q_2combs = list(combinations(np.arange(K_Q),2))

    best_obj = np.inf
    idxQ2P = []

    for pair in Q_2combs:
        
        # remove matrix
        old_idx = np.arange(K_Q)
        old_idx = np.delete(old_idx, [pair[0],pair[1]])
        
        # update matrix
        new_idx = Q.shape[1]-2
        Q_comb = np.hstack([Q[:,old_idx],np.expand_dims(Q[:,pair[0]]+Q[:,pair[1]],1)])

        opt_obj, idxQ2P_comb = align_ILP(P,Q_comb)
        idxQ2P_comb = np.array(idxQ2P_comb)

        if opt_obj<best_obj:
            best_obj = opt_obj
            idxQ2P = np.zeros(K_Q)

            idxQ2P[old_idx] = idxQ2P_comb[np.arange(len(old_idx))]
            idxQ2P[pair[0]] = idxQ2P_comb[new_idx]
            idxQ2P[pair[1]] = idxQ2P_comb[new_idx]
            idxQ2P = list(idxQ2P.astype(int))
            
    return best_obj, idxQ2P



def alignQ_wrtP(P: np.ndarray, Q: np.ndarray, idxQ2P: List[int], merge=True) -> Tuple[np.ndarray,list[int]]:
    assert P.shape[0] == Q.shape[0], "P and Q must have the same number of rows (individuals)."
    assert P.shape[1] <= Q.shape[1], "Q must have no fewer rows (clusters) than P."
    
    def default_diff_and_idx() -> Tuple[float, int]:
        return (float('inf'), -1)

    if merge:      
        aligned_Q = np.zeros_like(P)  
        new_pattern = idxQ2P
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,idxQ2P[q_idx]] += Q[:,q_idx]
    else:
        aligned_Q = np.zeros_like(Q)
        # find multiple Q columns that map to the same P column
        dups = np.unique([i for i in idxQ2P if idxQ2P.count(i)>1])
        extras = list()
        dups_min = defaultdict(default_diff_and_idx)

        new_pattern = [0 for _ in range(Q.shape[1])]
        for q_idx in range(Q.shape[1]):
            p_idx = idxQ2P[q_idx]
            if p_idx not in dups:
                new_pattern[q_idx] = p_idx
                aligned_Q[:,p_idx] = Q[:,q_idx]
            else:
                diff = np.linalg.norm(Q[:,q_idx]-P[:,p_idx])
                if dups_min[p_idx][0] > diff:
                    dups_min[p_idx] = (float(diff),q_idx) 
        # assign duplicated Q columns (except for the best-matched one) to new columns in aligned_Q
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
            aligned_Q[:,P.shape[1]+ie] = Q[:,e]

    validate_membership_matrix(aligned_Q)
    return aligned_Q, new_pattern


def align_within_k(Q_list: List[np.ndarray]):
    out_summary = []
    n_ids = len(Q_list)
    if n_ids < 2:
        logger.warning("Only one Q file is founded. No alignment within K clusters is performed.")
        return out_summary
    for i in range(n_ids-1):
        P = Q_list[i]
        for j in range(i+1,n_ids):
            Q = Q_list[j]
            opt_obj, idxQ2P = align_ILP(P, Q)
            # idxP2Q = reverse_alignment_same_K(idxQ2P)
            
            # calculate cost
            aligned_Q, _ = alignQ_wrtP(P,Q,idxQ2P,merge=True)
            cost = cost_membership(P, aligned_Q)
            
            # save output info
            out_summary.append([(i,j), cost, idxQ2P])

    return out_summary


def align_within_k_all_K(Q_names: List[str], K_range: List[int], K2IDs: dict,
                         qfile_dir: str,
                         output_dir: str, within_k_subdir: str = "alignment_withinK") -> Tuple[list, list]:
    """ Align Q files within each K value and save results.
        Args:
            Q_names: List of Q file base names (without extensions).
            K_range: List of unique K values.
            K2IDs: Dictionary mapping K values to lists of IDs.
            qfile_dir: Directory containing Q files.
            output_dir: Directory to save the output results.
            within_k_subdir: Subdirectory name for within-K alignment results.
        Returns:
            alignment_withinK_list: List of dictionaries mapping pairs of Q names to their alignment patterns.
            cost_withinK_list: List of dictionaries mapping pairs of Q names to their alignment costs.
    """
    within_k_dir = os.path.join(output_dir,within_k_subdir)
    if os.path.exists(within_k_dir) and os.listdir(within_k_dir):
        shutil.rmtree(within_k_dir)
        logger.info(f"Within-K output directory '{within_k_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(within_k_dir, exist_ok=True)
    logger.info(f"Created within-K output directory '{within_k_dir}'.")

    alignment_withinK_list = list()
    cost_withinK_list = list()
    for K in K_range:
        ids = K2IDs[K]
        # load Q files with the same K
        Q_subnames = [qname for qname in Q_names if qname.startswith(tuple(f"{i}_" for i in ids))]
        Q_files = [os.path.join(qfile_dir, f"{qname}.Q") for qname in Q_subnames]
        Q_list = load_withinK_qfiles(Q_files, K=K)  # Example K value, replace with actual logic to determine K
        disp_msg(f"Aligning {len(Q_list)} Q files within K={K}")
        alined_results = align_within_k(Q_list)
        # save results to file
        disp_msg(f"Saving alignment results within K={K}")
        write_aligned_within_k(alined_results, Q_subnames, os.path.join(within_k_dir,"K{}.txt".format(K)))
        # save results to dict
        alignment_withinK, cost_withinK = aligned_res_to_dicts(alined_results, Q_subnames)
        alignment_withinK_list.append(alignment_withinK)
        cost_withinK_list.append(cost_withinK)
        logger.info(f"{len(alignment_withinK)}({len(cost_withinK)}) pairwise results of alignment within K={K} completed.")

    return alignment_withinK_list, cost_withinK_list


def detect_modes_all_K(K_range, cost_withinK_list, Q_names, K2IDs, 
                       test_comm: bool = True, method: str = "louvain", res: float = 1.0, 
                       comm_min: float = 1e-4, comm_max: float = 1e-2) -> Tuple[list, list]:
    """ Detect modes for all K values based on the cost matrices.
        Args:
            K_range: List of unique K values.
            cost_withinK_list: List of dictionaries mapping pairs of Q names to their alignment costs.
            Q_names: List of Q file base names (without extensions).
            K2IDs: Dictionary mapping K values to lists of IDs.
            test_comm: Boolean indicating whether to test for community structure.
            method: Community detection method to use.
            res: Resolution parameter for the default Louvain community detection method.
            comm_min: Minimum threshold for the cost matrix.
            comm_max: Maximum threshold for the cost matrix.
        Returns:
            modes_all_K: List of detected modes for each K value.
            cost_matrices: List of corresponding cost matrices for each K value.
    """     
    modes_all_K = list()
    cost_matrices = list()
    for i_K, K in enumerate(K_range):
        disp_msg(f"Detecting modes for K={K}")

        cost_withinK = cost_withinK_list[i_K]
        # get Q names for this K
        Q_names_K = [qname for qname in Q_names if qname.startswith(tuple(f"{i}_" for i in K2IDs[K]))]
        assert len(Q_names_K) == len(K2IDs[K]), f"Number of replicate names ({len(Q_names_K)}) does not match number of IDs ({len(K2IDs[K])}) for K={K}."
        # construct cost matrix
        cost_mat = construct_cost_mat(cost_withinK, Q_names_K) 
        cost_matrices.append(cost_mat)

        # detect communities
        communities = detect_communities(cost_mat, test_comm = test_comm, method = method, res = res, 
                            min_threshold = comm_min, max_threshold = comm_max)
        modes = community_labels_to_modes(communities)
        modes_all_K.append(modes)
        logger.info(f"Detected {len(modes)} modes for K={K}.")
    
    return modes_all_K, cost_matrices


def cd_default(adj_mat: np.ndarray, res: float = 1.0, method: str = "louvain") -> List[int]:
    
    n_nodes = adj_mat.shape[0]
    G = nx.from_numpy_array(adj_mat)   
    G.remove_edges_from(nx.selfloop_edges(G))

    if method == "custom":
        logger.info(f"Using custom community detection method.")
        communities = cd_custom(adj_mat)
        return communities

    if method == "louvain":
        coms = algorithms.louvain(G, resolution=res)
        logger.info(f"Louvain detected {len(coms.communities)} communities with resolution {res}.")
    else:
        method_map = {
            "louvain": algorithms.louvain,
            "leiden": algorithms.leiden,
            "infomap": algorithms.infomap,
            "markov_clustering": algorithms.markov_clustering,
            "label_propagation": algorithms.label_propagation,
            "walktrap": algorithms.walktrap,
        }
        algr = method_map.get(method.lower())
        if algr is None:
            raise ValueError(f"Unknown method: {method}. Supported: {list(method_map.keys())}")
        coms = algr(G)
        logger.info(f"{method.capitalize()} detected {len(coms.communities)} communities.")
    
    communities = [-1] * n_nodes
    for cluster_idx, nodes in enumerate(coms.communities):
        for node in nodes:
            communities[node] = cluster_idx

    return communities


def detect_communities(cost_mat: np.ndarray, test_comm: bool = True, method: str = "louvain", res: float = 1.0, 
                       min_threshold: float = 1e-4, max_threshold: float = 1e-2) -> List[int]:
    """ Detect communities from a cost matrix.
        Args:
            cost_mat: A square cost matrix (2D NumPy array).
            threshold: A float value to threshold the cost matrix.
            norm: Boolean indicating whether to normalize the cost matrix.
        Returns:
            A binary adjacency matrix (2D NumPy array) where entries are 1 if the cost is below the threshold, else 0.
    """
    n_nodes = cost_mat.shape[0]
    if n_nodes == 0:
        raise ValueError("Cost matrix is empty.")
    if n_nodes == 1:
        logger.warning("Cost matrix has only one replicate. Returning a single singleton mode.")
        return [0]
    
    mask = ~np.eye(cost_mat.shape[0], dtype=bool) # mask for non-diagonal elements
    cost_min = np.min(cost_mat[mask]) 
    cost_max = np.max(cost_mat[mask])
    
    if n_nodes == 2:
        logger.warning("Cost matrix has two replicates. Returning modes based on thresholds.")
        c = cost_mat[0, 1]
        if c < min_threshold:
            logger.info(f"Cost {c} is smaller than threshold {min_threshold}. Returning a single mode.")
            return [0, 0]
        else:
            logger.info(f"Cost {c} is larger than threshold {min_threshold}. Returning two singleton modes.")
            return [0, 1]

    # check if all costs are small or large  
    if cost_max < min_threshold:
        logger.info(f"Largest cost {cost_max} is smaller than threshold {min_threshold}. Returning a single mode.")
        return [0]*n_nodes    
    if cost_min > max_threshold:
        logger.info(f"Smallest cost {cost_min} is larger than threshold {min_threshold}. Returning all singleton modes.")
        return list(range(n_nodes))
    
    # when n_nodes > 2, we run community detection
    adj_mat = cost_to_adj(cost_mat, norm=True)
    has_comm_struc = True
    if test_comm:
        # test for community structure
        has_comm_struc = test_comm_struc(adj_mat, alpha = 0.01)
    if has_comm_struc:
        logger.info("Detected community structure in the cost matrix. Running community detection.")
        communities = cd_default(adj_mat, method=method, res=res)
        if not communities:
            logger.warning("No communities detected. Returning singleton modes.")
            return list(range(n_nodes))
        # reorder community indices so that largest community is first, etc.
        counts = Counter(communities)
        # Sort indices by frequency (descending), then by first appearance
        order = sorted(counts.keys(),key=lambda x: (-counts[x], communities.index(x)))
        # Map old indices to new indices
        mapping = {old: new for new, old in enumerate(order)}
        # Apply mapping
        reordered_communities = [mapping[c] for c in communities]
        return reordered_communities
    else:
        logger.info("No significant community structure detected. Returning a single mode.")
        return [0]*n_nodes   


def find_repr_modes(modes: List[List[int]], cost_mat: np.ndarray):
    """ Find representative modes from community labels.
        Args:
            modes: A list of list of nodes in each mode.
            cost_mat: A square cost matrix (2D NumPy array).
        Returns:
            A list of representative modes.
    """
    repr_modes = []
    for mode_idx in range(len(modes)):
        all_indices = modes[mode_idx]
        assert isinstance(all_indices, list), "Expected modes to be a list of lists."
        if len(all_indices) == 1:
            repr_modes.append(all_indices[0])
        else:
            comm_cost_mat = cost_mat[np.array(all_indices), :][:, np.array(all_indices)]
            minC_idx = np.argmin(comm_cost_mat.sum(axis=0))
            repr_modes.append(all_indices[minC_idx])
            # other_rep_indices = [m for im, m in enumerate(all_indices) if im != minC_idx]
    return repr_modes


def extract_modes_and_stats(modes: List[List[int]], cost_mat: np.ndarray, 
                            Q_names: List[str], Q_list: List[np.ndarray], 
                            alignment_withinK: dict, label_prefix: str = "M") -> ModesDict:
    """ Extract modes and statistics from community labels.
        Args:
            modes: A list of list of nodes in each mode.
            cost_mat: A square cost matrix (2D NumPy array).
            alignment_withinK: A dictionary mapping pairs of Q names to their alignment indices.
            Q_names: A list of replicate names.
            Q_list: A list of Q matrices corresponding to the Q names.
            output_path: Directory to save the output results.
    """

    # extract modes and representative modes
    repr_modes = find_repr_modes(modes, cost_mat)
    assert len(modes) == len(repr_modes), "Number of modes and representative modes should match."
    
    repQ_modes = []
    avgQ_modes = []
    mode_labels = []
    mode_stats = []
    align_info = []
    for mode_idx in range(len(modes)):

        mode_label = "{}{}".format(label_prefix,mode_idx+1)
        mode_labels.append(mode_label)

        all_indices = modes[mode_idx]
        comm_cost_mat = cost_mat[np.array(all_indices),:][:,np.array(all_indices)]
        Gprime_mat = C2Gprime(comm_cost_mat)
        i_ref = repr_modes[mode_idx]
        repr_name = Q_names[i_ref]
        Q_mode = Q_list[i_ref]
        # if mode_idx == 0:
        #     major_mode = Q_mode
        #     major_mode_rep = repr_name
        # else:
        #     aligned_idxQ2P = alignment_withinK[(Q_names[i_ref],Q_names[r])]
        #     aligned_Q, _ = alignQ_wrtP(major_mode,Q,aligned_idxQ2P,merge=True)
        
        if len(all_indices)==1:
            mode_stats.append([mode_label, repr_name, 1, 0, 1])                    
            repQ_modes.append(Q_mode)
            avgQ_modes.append(Q_mode)
        else:
            # compute the average cost/similarity
            mode_stats.append([mode_label, repr_name, len(all_indices), avg_tril(comm_cost_mat), avg_tril(Gprime_mat)])
            # align this mode Q to prev mode Q
            repQ_modes.append(Q_mode)
            # align other Q's to the representative one
            # and get average Q over all reps in the mode
            Q_sum = np.zeros_like(Q_mode)
            Q_ref = Q_mode
            for r in all_indices:
                if r != i_ref:
                    Q = Q_list[r]
                    aligned_idxQ2P = alignment_withinK[(Q_names[repr_modes[mode_idx]],Q_names[r])]
                    # mode_alignment[Q_names[r]] = aligned_idxQ2P
                    align_info.append([mode_label, repr_name, Q_names[r], pattern_to_str(aligned_idxQ2P)])
                    aligned_Q, _ = alignQ_wrtP(Q_ref,Q,aligned_idxQ2P,merge=True)
                    Q_sum += aligned_Q
                else:
                    Q_sum += Q_ref
                    align_info.append([mode_label, repr_name, Q_names[i_ref], pattern_to_str(list(np.arange(Q_ref.shape[1])))])
                    # mode_alignment[Q_names[i_ref]] = np.arange(Q_ref.shape[1])
            avgQ_modes.append(Q_sum/len(all_indices))

    mode_stats = pd.DataFrame(mode_stats, columns=['Mode', 'Representative', 'Size', 'Cost', 'Performance'])
    align_info = pd.DataFrame(align_info, columns=['Mode', 'Representative', 'Replicate', 'Alignment'])
    
    return {"modes": modes, "repr_modes": repr_modes, \
            "mode_stats": mode_stats, "align_info": align_info,
            "repQ_modes": repQ_modes, "avgQ_modes": avgQ_modes}
    

def extract_modes_all_K(K_range: List[int], K2IDs: dict, Q_names: List[str], 
                        cost_matrices_list: List[np.ndarray], 
                        modes_all_K_list: List[list],
                        alignment_withinK_list: List[dict],
                        processed_input_dir: str, output_dir: str) -> Tuple[list, pd.DataFrame]:
    cd_res = []
    for i_K, K in enumerate(K_range):
        Q_names_K = [qname for qname in Q_names if qname.startswith(tuple(f"{i}_" for i in K2IDs[K]))]
        cost_mat = cost_matrices_list[i_K]
        modes = modes_all_K_list[i_K]
        Q_files_K = [os.path.join(processed_input_dir, f"{name}.Q") for name in Q_names_K]
        Q_list_K = load_withinK_qfiles(Q_files_K, K=K)
        cd_res_K = extract_modes_and_stats(modes, cost_mat, Q_names_K, 
                                           Q_list_K, alignment_withinK_list[i_K], 
                                           label_prefix="K{}M".format(K))
        cd_res.append(cd_res_K)

    # save mode detection results
    if os.path.exists(output_dir) and os.listdir(output_dir):
        shutil.rmtree(output_dir)
        logger.info(f"Mode output directory '{output_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(output_dir, exist_ok=True)
    
    # save results to files
    mode_stats = list()
    for cd_res_K in cd_res:
        mode_stats.append(cd_res_K['mode_stats'])
    mode_stats = pd.concat(mode_stats, ignore_index=True)
    mode_stats.to_csv(os.path.join(output_dir, "mode_stats.txt"), index=False)
    align_info = list()
    for cd_res_K in cd_res:
        align_info.append(cd_res_K['align_info'])
    align_info = pd.concat(align_info, ignore_index=True)
    align_info.to_csv(os.path.join(output_dir, "mode_alignment.txt"), index=False)
    avg_stat = list()

    for i_K, K in enumerate(K_range):
        cd_res_K = cd_res[i_K]
        # save mode memberships
        for mode_idx, mode_label in enumerate(cd_res_K['mode_stats']['Mode'].values):
            np.savetxt(os.path.join(output_dir,"{}_rep.Q".format(mode_label)), cd_res_K['repQ_modes'][mode_idx], delimiter=' ')
            np.savetxt(os.path.join(output_dir,"{}_avg.Q".format(mode_label)), cd_res_K['avgQ_modes'][mode_idx], delimiter=' ')    
        # compute average statistics
        avg_stat_K = compute_mode_avg_stats(cd_res_K)
        avg_stat_K['K'] = K
        avg_stat_K = avg_stat_K[['K', 'Size', 'NS-Size', 'Cost', 'Performance']]
        avg_stat.append(avg_stat_K)
    # save average statistics
    avg_stat = pd.concat(avg_stat, ignore_index=True)
    avg_stat.to_csv(os.path.join(output_dir, "mode_average_stats.txt"), index=False) 

    return cd_res, avg_stat


def align_across_k(K_range: List[int], Q_list_list: List[List[np.ndarray]], Q_labels_list: List[List[str]], merge: bool=True) -> Tuple[dict, dict, pd.DataFrame, pd.DataFrame]:
    """ Align Q files across different K values.
        Args:
            K_range: List of unique K values.
            Q_list_list: List of lists of Q matrices for each K value.
            Q_labels_list: List of lists of labels for each Q matrix.
            merge: Boolean indicating whether to merge clusters in the alignment.
        Returns:
            alignment_acrossK: Dictionary mapping pairs of Q names to their alignment indices.
            cost_acrossK: Dictionary mapping pairs of Q names to their alignment costs.
            best_acrossK_out_summary: DataFrame summarizing the best alignments across K values.
            major_acrossK_out_summary: DataFrame summarizing the major alignments across K values.
    """
    # Check input consistency
    assert len(K_range) == len(Q_list_list), "K_range does not match the length of Q_list_list."
    assert len(K_range) == len(Q_labels_list), "K_range does not match the length of Q_labels_list."
    
    K_range_sorted = sorted(K_range,reverse=True) # reverse order: from large K to small K
    K_comb = list([(K_range_sorted[i],K_range_sorted[i+1]) for i in range(len(K_range_sorted)-1)])
    K_comb.extend([(K_range_sorted[i],K_range_sorted[i]) for i in range(len(K_range_sorted))])

    alignment_acrossK = dict()
    cost_acrossK = dict()
    best_alignments = dict()

    for K1,K2 in K_comb:
        disp_msg(f"Aligning K={K1} to K={K2} modes")
        
        i_K1 = K_range.index(K1)
        i_K2 = K_range.index(K2)

        best_alignment_idx = 0
        best_alignment_obj = np.inf

        # all pairs of Q modes across K1 and K2
        rijs = list(product(range(len(Q_list_list[i_K1])),range(len(Q_list_list[i_K2]))))
        
        for i,(ri,rj) in enumerate(rijs):
            
            Q1 = Q_list_list[i_K1][ri] # Q1 is the one with more clusters
            Q2 = Q_list_list[i_K2][rj]

            if (K1==K2 and rj==ri):
                # identical alignment (dummy)
                opt_obj = 0
                idxQ2P = np.arange(K1)
                idxQ2P = list(idxQ2P.astype(int))
            elif merge and (K1-K2)==1:
                opt_obj, idxQ2P = align_ILP_diff1(Q2, Q1)
            else:
                opt_obj, idxQ2P = align_ILP(Q2, Q1)

            # check if this alignment is the best one between K1 and K2
            if best_alignment_obj>opt_obj:
                best_alignment_obj = opt_obj
                best_alignment_idx = i
            # calculate cost
            aligned_Q, _ = alignQ_wrtP(Q2,Q1,idxQ2P,merge=True)
            cost = cost_membership(Q2,aligned_Q)
            lb2 = Q_labels_list[i_K2][rj]
            lb1 = Q_labels_list[i_K1][ri]
            pair_label = "{}-{}".format(lb2, lb1)
            
            cost_acrossK[pair_label] = cost
            alignment_acrossK[pair_label] = idxQ2P    
        best_alignments[(K1,K2)] = rijs[best_alignment_idx]

    # best alignments
    best_acrossK_out_summary = []
    for ii in range(len(K_range_sorted)-1):
        K1 = K_range_sorted[ii]
        K2 = K_range_sorted[ii+1]
        i_K1 = K_range.index(K1)
        i_K2 = K_range.index(K2)
        ri,rj = best_alignments[(K1,K2)]
        pair_label = "{}-{}".format(Q_labels_list[i_K2][rj], Q_labels_list[i_K1][ri])
        Q1 = Q_list_list[i_K1][ri] 
        Q2 = Q_list_list[i_K2][rj]
        cost_pair = cost_acrossK[pair_label]
        cost_sep = cost_membership_sep(Q2,Q1,alignment_acrossK[pair_label])
        best_acrossK_out_summary.append([pair_label,cost_pair,C2Gprime(np.array([cost_pair]))[0],cost_sep,C2Gprime(np.array([cost_sep]))[0]])
    best_acrossK_out_summary = pd.DataFrame(best_acrossK_out_summary, columns=['Best Pair', 'Cost', 'Performance', 'Separate-Cluster Cost', 'Separate-Cluster Performance'])    

    # major alignments
    major_acrossK_out_summary = []
    for ii in range(len(K_range_sorted)-1):
        K1 = K_range_sorted[ii]
        K2 = K_range_sorted[ii+1]
        i_K1 = K_range.index(K1)
        i_K2 = K_range.index(K2)
        ri,rj = 0,0
        pair_label = "{}-{}".format(Q_labels_list[i_K2][rj],Q_labels_list[i_K1][ri])
        Q1 = Q_list_list[i_K1][ri] 
        Q2 = Q_list_list[i_K2][rj]
        cost_pair = cost_acrossK[pair_label]
        cost_sep = cost_membership_sep(Q2,Q1,alignment_acrossK[pair_label])
        major_acrossK_out_summary.append([pair_label,cost_pair,C2Gprime(np.array([cost_pair]))[0],cost_sep,C2Gprime(np.array([cost_sep]))[0]])
    major_acrossK_out_summary = pd.DataFrame(major_acrossK_out_summary, columns=['Major Pair', 'Cost', 'Performance', 'Separate-Cluster Cost', 'Separate-Cluster Performance'])    

    return alignment_acrossK, cost_acrossK, best_acrossK_out_summary, major_acrossK_out_summary


def write_alignment_across_k(alignment_acrossK: dict, cost_acrossK: dict, output_file: str):
    """ Write the alignment across K results to a file.
        Args:
            alignment_acrossK: Dictionary mapping pairs of Q names to their alignment indices.
            cost_acrossK: Dictionary mapping pairs of Q names to their alignment costs.
            out_summary: List of lists containing pairwise results of alignment across K.
            output_path: Path to the output file.
    """
    with open(output_file, "w") as f:
        f.write("Mode1-Mode2,Cost,Alignment\n")
        for pair_label in alignment_acrossK.keys():
            idxQ2P = alignment_acrossK[pair_label]
            cost = cost_acrossK[pair_label]
            f.write("{},{},{}\n".format(pair_label, cost, pattern_to_str(idxQ2P)))

# Reorder clusters in Q matrices according to the alignment pattern
def reorderQ_within_k(Q_list: list[np.ndarray], Q_names: list[str], alignment_acrossK: dict) -> list[np.ndarray]:
    """
    Reorder clusters in Q matrices within a single K value.
    
    Args:
        Q_list: List of matrices loaded from Q files.
        Q_names: List of names corresponding to the Q matrices.
        alignment_acrossK: Dictionary containing alignment patterns.
        
    Returns:
        List of aligned Q matrices.
    """
    if not Q_list:
        raise ValueError("Q_list is empty. Cannot align Q matrices.")
    
    aligned_Qs = []
    lb1 = Q_names[0]  # Assuming all Q matrices have the same label for within-K alignment
    Q = Q_list[0]
    aligned_Qs.append(Q)
    
    for i in range(1,len(Q_list)):
        lb2 = Q_names[i]
        # retrieve alignment pattern
        alignment = alignment_acrossK["{}-{}".format(lb1, lb2)]
        Q = Q_list[i]
        aligned_Q, _ = alignQ_wrtP(Q,Q,alignment,merge=True)
        aligned_Qs.append(aligned_Q)
    return aligned_Qs

def reorderQ_across_k(K_range: list[int], Q_modes_list: list[list[np.ndarray]], mode_names_list: list[list[str]], 
                      alignment_acrossK: dict, anchor_pairs: list[str]) -> Tuple[dict, dict]:
    
    anchor_pairs_rev = anchor_pairs[::-1]
    all_modes_alignment = {lb:[] for i_K, K in enumerate(K_range) for lb in mode_names_list[i_K]}
    aligned_Qs_allK = dict()
    base_patterns = dict()
    i_K = 0
    K = K_range[i_K]
    mode_names = mode_names_list[i_K]
    if len(anchor_pairs)>0:
        m1 = anchor_pairs_rev[0].split("-")[0]
        m_m1 = mode_names.index(m1)
        base_Q = Q_modes_list[i_K][m_m1]
        base_patterns[m1] = [i for i in range(K)]
        all_modes_alignment[m1] = [i for i in range(base_Q.shape[1])]
    else:
        m_m1 = 0
        m1 = mode_names[m_m1]
        base_Q = Q_modes_list[i_K][0]
        all_modes_alignment[m1] = list(range(base_Q.shape[1]))
    aligned_Qs_allK[m1] = base_Q
    
    modes_indices = range(len(Q_modes_list[i_K]))
    if len(modes_indices)>1:
        for m in modes_indices:
            if m!=m_m1:
                m2 = mode_names_list[i_K][m]
                # retrieve alignment pattern
                ali_pat = alignment_acrossK["{}-{}".format(m1,m2)]
                Q = Q_modes_list[i_K][m]
                aligned_Q = np.zeros_like(Q)
                for q_idx in range(Q.shape[1]):
                    aligned_Q[:,ali_pat[q_idx]] += Q[:,q_idx]
                base_patterns[m2] = ali_pat
                all_modes_alignment[m2] = ali_pat
                aligned_Qs_allK[m2] = aligned_Q

    for i_pair, pair in enumerate(anchor_pairs_rev):
        m1 = pair.split("-")[0]
        m2 = pair.split("-")[1]
        i_K = i_pair + 1
        K = K_range[i_K]
        modes_indices = range(len(Q_modes_list[i_K]))
        m_m2 = mode_names_list[i_K].index(m2)
        Q = Q_modes_list[i_K][m_m2]
        m1_K = K_range[i_pair]
        P = Q_modes_list[i_pair][mode_names_list[i_pair].index(m1)]
        pattern = alignment_acrossK[pair]
        aligned_Q, new_pattern = alignQ_wrtP(P,Q,pattern,merge=False)
        pat = [base_patterns[m1][i] if i<m1_K else i for i in new_pattern]
        aligned_Q = np.zeros_like(aligned_Q)
        for q_idx in range(Q.shape[1]):
            aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
        base_patterns[m2] = pat
        all_modes_alignment[m2] = pat
        aligned_Qs_allK[m2] = aligned_Q
        
        if len(modes_indices)>1:
            for m in modes_indices:
                if m!=m_m2:
                    m3 = mode_names_list[i_K][m]
                    ali_pat = alignment_acrossK["{}-{}".format(m2,m3)]
                    Q = Q_modes_list[i_K][m]
                    P = Q_modes_list[i_K][m_m2]
                    aligned_Q, new_pattern = alignQ_wrtP(P,Q,ali_pat,merge=False)
                    
                    pat = [base_patterns[m2][i] for i in new_pattern ]
                    aligned_Q = np.zeros_like(aligned_Q)
                    for q_idx in range(Q.shape[1]):
                        aligned_Q[:,pat[q_idx]] += Q[:,q_idx]
                    base_patterns[m3] = pat
                    all_modes_alignment[m3] = pat
                    aligned_Qs_allK[m3] = aligned_Q
    
    for m in all_modes_alignment.keys():
        orig =  all_modes_alignment[m]
        new = [orig.index(ii) for ii in range(len(orig))]
        all_modes_alignment[m] = new

    return aligned_Qs_allK, all_modes_alignment
