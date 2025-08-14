import os
import shutil
import numpy as np
import argparse
import importlib.resources
from typing import List, Any, Tuple, Optional #, Optional, Tuple
import matplotlib.colors as mcolors

import logging
logger = logging.getLogger(__name__)

def str2bool(v):
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError("Boolean value expected.")
    

def disp_params(args: argparse.Namespace, title: str = ""):
    logger.info(f" {title} ".center(50, '='))
    logger.info("--------------------------------------------------") 
    logger.info("Parameters:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")
    logger.info("--------------------------------------------------") 

def disp_msg(msg: str, level: str = "info"):
    s = f">>> {msg} ".ljust(50, '.')
    if level == "info":
        logger.info(s)
    elif level == "warning":
        logger.warning(s)
    elif level == "error":
        logger.error(s)
    else:
        logger.info(s)

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


def load_matrix(file_path: str, delimiter: str = " ", skip_rows: int = 0) -> np.ndarray:
    """
    Load a numeric matrix from a file, skipping missing data and specified rows.
    
    Args:
        file_path: Path to the input file.
        delimiter: Delimiter used in the file (default is space).
        skip_rows: Number of rows to skip at the top of the file.
        
    Returns:
        A NumPy array containing the loaded matrix.
    """
    matrix = np.loadtxt(file_path, delimiter=delimiter, skiprows=skip_rows, dtype=float)
    return matrix


def cost_membership(P: np.ndarray, Q: np.ndarray, ord: Any = 'fro') -> float:
    """ Compute the cost (squared Euclidean distance) between two membership matrices P and Q.
        Args:
            P: Membership matrix of shape (N_ind, K_P).
            Q: Membership matrix of shape (N_ind, K_Q).
            ord: Order of the norm to use as in numpy.linalg.norm (default is 'fro' for Frobenius norm).
        Returns:
            The cost value as a float, computed using the norm with designated order.
    """
    assert P.shape == Q.shape, "P and Q must have the same size."
    N_ind = P.shape[0]
    return float(np.linalg.norm(P-Q, ord=ord)**2/(2*N_ind))


def cost_membership_sep(P: np.ndarray, Q: np.ndarray, idxQ2P: List[int]) -> float:
    """ Compute the average of column-wise cost between two membership matrices P and Q
        Args:
            P: Membership matrix of shape (N_ind, K_P).
            Q: Membership matrix of shape (N_ind, K_Q).
            idxQ2P: List of indices mapping Q to P.
        Returns:
            The cost value as a float, computed using the norm with designated order.
    """
    N_ind = P.shape[0]
    cost_col = list()
    for p_idx in range(P.shape[1]):
        cost_col.append(list())
    for q_idx in range(Q.shape[1]):
        p_idx = idxQ2P[q_idx]
        cost_col[p_idx].append([np.sum((P[:,p_idx]-Q[:,q_idx])**2)])
    cost = 0    
    for p_idx in range(P.shape[1]):
        if len(cost_col[p_idx])==1:
            cost += np.mean(cost_col[p_idx])
    cost /= 2*N_ind
    return cost


def pattern_to_str(patt: List[int]) -> str:
    """ Convert a pattern (of index matching) list to a string representation.
        Args:
            perm: List of integers representing the index matching.
        Returns:
            A string representation, e.g., "1 2 3".
    """
    return " ".join(str(i + 1) for i in patt)


def str_to_pattern(s: str) -> List[int]:
    """ Convert a string representation of a pattern to a list of integers.
        Args:
            s: A string representation, e.g., "1 2 3".
        Returns:
            A list of integers representing the index matching.
    """
    return [int(i) - 1 for i in s.split()]

def reverse_alignment_same_K(idxQ2P: List[int]) -> List[int]:
    """ Reverse the alignment mapping from Q to P (when both have the same K).
        Example:
        If idxQ2P = [1, 2, 0] means Q[0] -> P[1], Q[1] -> P[2], Q[2] -> P[0]
        then idxP2Q should be [2, 0, 1].
    Args:
        idxQ2P: List of indices mapping Q to P.
    Returns:
        List of indices mapping P to Q.
    """
    idxP2Q = np.zeros(len(idxQ2P)).astype(int)
    for i,idx in enumerate(idxQ2P):
        idxP2Q[idx] = i
    return list(idxP2Q)


def parse_strings(strlist: str = "", strs: List[str]=[], remove_dup: bool=True) -> List[str]:
    """ Parse input strings from command-line arguments or a string list. 
    Args:
        strlist: Path to a plain text file containing strings (one per line) (optional).
        strs: A list of strings passed directly (space separated) (optional).
    Returns:
        A list of input file paths.
    """
    strings = []

    # Get strings from strlist if provided
    if strlist:
        if not os.path.isfile(strlist):
            raise FileNotFoundError(f"File list not found: {strlist}")
        with open(strlist) as f:
            strings_from_file = [line.strip() for line in f if line.strip()]
        strings.extend(strings_from_file)

    # Append any strings passed directly
    if strs:
        strings.extend(strs)

    if remove_dup:
        # Remove duplicates if needed while preserving order
        strings = list(dict.fromkeys(strings))  
    return strings



def C2Gprime(C: np.ndarray) -> np.ndarray:       
    """ Convert cost to Gprime.
    """
    return 1-np.sqrt(C)


def avg_tril(A: np.ndarray) -> float:
    """ Compute the average of the lower-triangular part of a square matrix (excluding diagonal
        Args:
            A: A square matrix (2D NumPy array).
        Returns:
            The average of the lower-triangular part of the matrix.
    """
    # assert A is symmetric
    assert A.shape[0] == A.shape[1], "Input matrix must be square."
    return np.sum(np.tril(A, k=-1))/(A.shape[0]*(A.shape[0]-1)/2)


def get_modes_all_K(K_range: List[int], cd_res: List[dict]) -> Tuple[List[List[str]], List[List[np.ndarray]], List[List[np.ndarray]]]:
    """ Extract mode names, representative modes, and average modes for all K values.
        Args:
            K_range: List of K values.
            cd_res: List of community detection results for each K.
        Returns:
            mode_names_list: List of lists of mode names for each K.
            Q_rep_modes_list: List of lists of representative mode matrices for each K.
            Q_avg_modes_list: List of lists of average mode matrices for each K.
    """
    # get modes for alignment across K
    mode_names_list = []
    Q_rep_modes_list = []
    Q_avg_modes_list = []
    for i_K, K in enumerate(K_range):
        mode_names = [mode_label for mode_label in cd_res[i_K]['mode_stats']['Mode'].values]
        mode_names_list.append(mode_names)
        Q_rep_modes = [cd_res[i_K]['repQ_modes'][mode_idx] for mode_idx in range(len(mode_names))]
        Q_avg_modes = [cd_res[i_K]['avgQ_modes'][mode_idx] for mode_idx in range(len(mode_names))]
        Q_rep_modes_list.append(Q_rep_modes)
        Q_avg_modes_list.append(Q_avg_modes)
    
    return mode_names_list, Q_rep_modes_list, Q_avg_modes_list


def write_reordered_across_k(aligned_Qs_allK: dict, all_modes_alignment: dict, output_dir: str, suffix: str):
    """ Write the aligned Q matrices across K values to files.
        Args:
            aligned_Qs_allK: Dictionary of aligned Q matrices for each mode.
            all_modes_alignment: Dictionary of alignment patterns for each mode.
            output_dir: Directory to save the aligned Q matrices.
    """
    assert len(aligned_Qs_allK)==len(all_modes_alignment), "Number of aligned Q matrices does not match number of alignment patterns."
    if os.path.exists(output_dir) and os.listdir(output_dir):
        shutil.rmtree(output_dir)
        logger.info(f"Aligned membership matrices output directory '{output_dir}' already exists and is not empty. Removed existing directory.")
    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, f"all_modes_alignment_{suffix}.txt"), "w") as f:
        for mode_label, aligned_Q in aligned_Qs_allK.items():
            np.savetxt(os.path.join(output_dir, f"{mode_label}_{suffix}.Q"), aligned_Q, delimiter=' ')
            alignment_pattern = all_modes_alignment[mode_label]
            f.write(f"{mode_label}:{pattern_to_str(alignment_pattern)}\n")

def get_mode_sizes(cd_res: list) -> dict:
    """ Get sizes of modes for each K value.
        Args:
            cd_res: List of dictionaries containing mode detection results for each K value.
            K_range: List of unique K values.
        Returns:
            Dictionary containing mode sizes for each mode.
    """
    mode_sizes = dict()
    for i_K, cd_res_K in enumerate(cd_res):
        mode_stats = cd_res_K['mode_stats']
        for _, row in mode_stats.iterrows():
            mode_sizes[row['Mode']] = row['Size']
    
    return mode_sizes


def unnest_list(list_of_lists: List[List[Any]]) -> List[Any]:
    """ Flatten a list of lists into a single list.
        Args:
            list_of_lists: A list containing sublists.
        Returns:
            A flattened list containing all elements from the sublists.
    """
    return [item for sublist in list_of_lists for item in sublist]


def labels_are_grouped(labels: list, uniq_labels: list):
    for lb in uniq_labels:
        positions = [i for i, val in enumerate(labels) if val == lb]
        if positions and (max(positions) - min(positions) + 1 != len(positions)):
            return False
    return True

def get_uniq_lb_sep(labels: list) -> Tuple[list, list, list]:
    """ Get unique labels and their separation indices.
        Args:
            labels: List of labels.
        Returns:
            uniq_lbs: Unique labels.
            uniq_lbs_indices: Indices (middle) of unique labels.
            uniq_lbs_sep_idx: Separation indices for unique labels.
    """
    uniq_lbs = list(np.unique(labels))
    # assert that items with same lbs are together:
    assert labels_are_grouped(labels, uniq_lbs), "Labels are not grouped together."
    uniq_lbs_indices = list()
    uniq_lbs_sep_idx = list()
    uniq_lbs_sep_idx.append(0)
    for lb in uniq_lbs:
        indices = [i for i, val in enumerate(labels) if val == lb]
        mean_idx = np.mean(indices)
        uniq_lbs_indices.append(mean_idx)
        uniq_lbs_sep_idx.append(indices[-1])
    return uniq_lbs, uniq_lbs_indices, uniq_lbs_sep_idx