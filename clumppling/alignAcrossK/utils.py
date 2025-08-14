import numpy as np
from typing import Tuple
from clumppling.utils import load_matrix
import logging  
logger = logging.getLogger(__name__)

def load_any_qfiles(qfiles: list[str]) -> list:
    """
    Load Q files (any K). 
    Args:
        qfiles: List of Q file paths. 
    Returns:
        List of matrices loaded from the Q files.
    """
    if not qfiles:
        raise ValueError("No Q files provided for loading.")
 
    Q_list = []
    for file_path in qfiles:
        try:
            matrix = load_matrix(file_path)
            Q_list.append(matrix)
        except Exception as e:
            logger.error(f"ERROR loading {file_path}: {e}")
    
    logger.info(f"Successfully loaded {len(Q_list)} Q files.")

    return Q_list


def extract_K_range_from_Qs(Q_list: list) -> list:
    """
    Extract unique K values from a list of Q matrices.
    
    Args:
        Q_list: List of matrices loaded from Q files.
        
    Returns:
        List of unique K values.
    """
    if not Q_list:
        raise ValueError("Q_list is empty. Cannot extract K range.")
    
    K_range = sorted(set(matrix.shape[1] for matrix in Q_list))
    logger.info(f"Extracted K range: {K_range}")
    
    return K_range


def separate_Qs_by_K(K_range: list[int], Q_list: list[np.ndarray], name_list: list[str]) -> Tuple[list,list]:
    """
    Separate Q matrices by their K values.
    
    Args:
        Q_list: List of matrices loaded from Q files.
        K_range: List of unique K values.
        
    Returns:
        Dictionary mapping each K value to a list of corresponding Q matrices.
    """
    Qs_by_K = {K: [] for K in K_range}
    names_by_K = {K: [] for K in K_range}
    
    for i in range(len(Q_list)):
        matrix = Q_list[i]
        K = matrix.shape[1]
        if K in Qs_by_K:
            Qs_by_K[K].append(matrix)
            names_by_K[K].append(name_list[i])
        else:
            logger.warning(f"Matrix with unsupported K value {K} found. Skipping.")
    
    # dict to list conversion
    Qs_by_K = [Qs_by_K[K] for K in K_range]
    names_by_K = [names_by_K[K] for K in K_range]
    return Qs_by_K, names_by_K