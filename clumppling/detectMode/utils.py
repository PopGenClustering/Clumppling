import os
import shutil
import numpy as np
import pandas as pd
from typing import List, Optional, TypedDict
from TracyWidom import TracyWidom

# from clumppling.core import aligned_res_to_dicts
# from clumppling.alignWithinK import load_aligned_within_k
import logging  
logger = logging.getLogger(__name__)

def construct_cost_mat(cost_withinK: dict, Q_names: List[str]) -> np.ndarray:
    cost_mat = np.zeros((len(Q_names), len(Q_names)))
    for key, value in cost_withinK.items():
        id1, id2 = key
        i, j = Q_names.index(id1), Q_names.index(id2)
        cost_mat[i, j] = value
        cost_mat[j, i] = value
    # check cost matrix validity
    if not np.all(np.isfinite(cost_mat)):
        raise ValueError("Cost matrix contains non-finite values (NaN or Inf).")
    if np.any(cost_mat < 0):
        raise ValueError("Cost matrix contains negative values.")
    if np.any(cost_mat > 1):
        raise ValueError("Cost matrix contains values above 1.")
    if not np.all(np.diag(cost_mat) == 0):
        raise ValueError("Cost matrix diagonal is not zero.")
    return cost_mat


def standardize_matrix(W: np.ndarray) -> np.ndarray:
    """ Standardize the input matrix W by subtracting the mean and dividing by the standard deviation.
        Args:
            W: A square matrix (2D NumPy array).
        Returns:
            A standardized matrix (2D NumPy array) where each off-diagonal element is standardized.
    """
    non_diag_idx = np.where(~np.eye(W.shape[0],dtype=bool))
    non_diag_w = W[non_diag_idx]
    W_standardized = np.zeros(W.shape)
    W_standardized[non_diag_idx] = (non_diag_w-non_diag_w.mean())/non_diag_w.std()
    return W_standardized

def normalize_matrix(W: np.ndarray) -> np.ndarray:
    """ Normalize the input matrix W by dividing each element by the square root of the number of rows.
        Args:
            W: A square matrix (2D NumPy array).
        Returns:
            A normalized matrix (2D NumPy array) where each element is divided by the square root of the number of rows.
    """
    return W/np.sqrt(W.shape[0])

def exponentiate_matrix(W: np.ndarray, t: float) -> np.ndarray:
    """ Exponentiate non-diagonal entries of the input matrix W with scale t.
        Args:
            W: A square matrix (2D NumPy array).
            t: A float value to exponentiate the off-diagonal elements.
        Returns:
            An exponentiated matrix (2D NumPy array) where each off-diagonal element is raised to the power of t.
    """
    non_diag_idx = np.where(~np.eye(W.shape[0],dtype=bool))
    W_exp = np.zeros(W.shape)
    W_exp[non_diag_idx] = np.exp(W[non_diag_idx]*t)
    return W_exp


def test_comm_struc(W: np.ndarray, alpha: float = 0.01):
    """ Test if the input matrix W has community structure using Tracy-Widom test.
        Args:
            W: A square matrix (2D NumPy array).
            alpha: Significance level for the test.
        Returns:
            True if the matrix has community structure, False otherwise.
    """
    # check if W is square
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError("Input matrix W must be square (2D NumPy array).")
    
    # standardization mapping
    W_standardized = standardize_matrix(W)
    T = normalize_matrix(W_standardized)
    
    W_exp = exponentiate_matrix(W, t = 1/2)
    Te = normalize_matrix(standardize_matrix(W_exp))
    
    # Tracy Widom CI
    x = np.linspace(-10, 10, 1001)
    tw = TracyWidom(beta=1)  # allowed beta values are 1, 2, and 4
    cdf = tw.cdf(x)
    
    CI_max = x[np.where(cdf>(1-alpha/4))[0][0]]
    CI_min = x[np.where(cdf<alpha/4)[0][-1]]
    CI_p_max, CI_p_min = CI_max, CI_min
    
    # test
    eigval, _ = np.linalg.eig(T)
    eig_T_max, eig_T_min = eigval.max(), eigval.min()
    eigval, _ = np.linalg.eig(Te)
    eig_Te_max, eig_Te_min = eigval.max(), eigval.min()
    
    s = int(eig_T_max<CI_max) + int(eig_T_min>CI_min) + int(eig_Te_max<CI_p_max) + int(eig_Te_min>CI_p_min)
    has_comm_struc = s!=4
    
    return has_comm_struc


def cost_to_adj(cost_mat: np.ndarray, norm: bool) -> np.ndarray:
    """ Convert cost matrix to adjacency matrix using a threshold.
        Args:
            cost_mat: A square cost matrix (2D NumPy array).
            threshold: A float value to threshold for small cost.
            norm: Boolean indicating whether to normalize the cost matrix.
        Returns:
            An adjacency matrix (2D NumPy array) where entries are 1 if the cost is 0.
    """
    # Note: all cost_mat entries should be in [0,1].
    if norm:
        mask = ~np.eye(cost_mat.shape[0], dtype=bool) # mask for non-diagonal elements
        cost_min = np.min(cost_mat[mask]) 
        cost_max = np.max(cost_mat[mask])
        adj_mat = 1-(cost_mat-cost_min)/(cost_max-cost_min)
    else:
        adj_mat = 1-cost_mat
    np.fill_diagonal(adj_mat, 1)  # ensure self-loops
    return adj_mat


def community_labels_to_modes(communities: List[int]) -> List[List[int]]:
    """ Convert community labels to modes.
        Args:
            communities: A list of community labels for each node.
        Returns:
            A list of modes, where each mode corresponds to list of node indices.
    """
    if not communities:
        return []
    modes = []
    assert max(communities)==len(np.unique(communities)) - 1, "Community labels should be zero-indexed and consecutive."
    for i in range(max(communities) + 1):
        mode = []
        for j, c in enumerate(communities):
            if c == i:
                mode.append(j)
        modes.append(mode)
    return modes


class ModesDict(TypedDict):
    modes: List[List[int]]
    repr_modes: List[int]
    mode_stats: pd.DataFrame
    align_info: pd.DataFrame
    repQ_modes: List[np.ndarray]
    avgQ_modes: List[np.ndarray]

def compute_mode_avg_stats(res: ModesDict) -> pd.DataFrame:
    """ Compute average statistics of modes.
        Args:
            res: Dictionary containing modes, representative modes, mode stats, alignment info, repQ modes, and avgQ modes.
        Returns:
            avg_stat: DataFrame containing average statistics of modes.
    """
    total_s = sum([len(r) for r in res['modes']])
    non_singleton_cost = 0
    non_singleton_perf = 0
    non_singleton_s = 0
    non_singleton_pairwise_cnt = 0
    for mode_idx in range(len(res['mode_stats'])):
        s = res['mode_stats'].iloc[mode_idx]["Size"]
        if s > 1:
            pairwise_cnt = s * (s - 1) // 2
            non_singleton_cost += res['mode_stats'].iloc[mode_idx]['Cost'] * pairwise_cnt
            non_singleton_perf += res['mode_stats'].iloc[mode_idx]['Performance'] * pairwise_cnt
            non_singleton_s += s
            non_singleton_pairwise_cnt += pairwise_cnt
    if non_singleton_s > 0:
        avg_stat = [total_s, non_singleton_s, 
                    non_singleton_cost / non_singleton_pairwise_cnt, 
                    non_singleton_perf / non_singleton_pairwise_cnt]
    else:
        avg_stat = [total_s, non_singleton_s, 0, 1]
    return pd.DataFrame([avg_stat], columns=['Size', 'NS-Size', 'Cost', 'Performance'])


def write_modes_to_file(res: ModesDict, output_dir: str, compute_avg: bool=True) -> Optional[pd.DataFrame]:
    """
    Write the results of mode detection to files in the specified output directory.
    
    Args:
        res: Dictionary containing modes, representative modes, mode stats, alignment info, repQ modes, and avgQ modes.
        output_dir: Directory to save the output result files.
        compute_avg: Whether to compute and save average statistics of modes (default: True).
    Returns:
        avg_stat: DataFrame containing average statistics of modes (optional).
    """

    # prepare output directory
    if os.path.exists(output_dir) and os.listdir(output_dir):
        shutil.rmtree(output_dir)
        logger.info(f"Mode output directory '{output_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(output_dir, exist_ok=True)
    
    # save results to files
    res['mode_stats'].to_csv(os.path.join(output_dir, "mode_stats.txt"), index=False)
    res['align_info'].to_csv(os.path.join(output_dir, "mode_alignments.txt"), index=False)
    for mode_idx, mode_label in enumerate(res['mode_stats']['Mode'].values):
        np.savetxt(os.path.join(output_dir,"{}_rep.Q".format(mode_label)), res['repQ_modes'][mode_idx], delimiter=' ')
        np.savetxt(os.path.join(output_dir,"{}_avg.Q".format(mode_label)), res['avgQ_modes'][mode_idx], delimiter=' ')    

    if compute_avg:
        logger.info(f"Compute average statistics over all modes.")
        # compute and save average statistics
        avg_stat = compute_mode_avg_stats(res)
        avg_stat.to_csv(os.path.join(output_dir, "mode_average_stats.txt"), index=False) 
        return avg_stat
    else:
        return None

