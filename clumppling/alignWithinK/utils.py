import numpy as np
from typing import List, Any, Tuple
from clumppling.utils import load_matrix, str_to_pattern, pattern_to_str, reverse_alignment_same_K
import logging  
logger = logging.getLogger(__name__)

def load_withinK_qfiles(qfiles: list[str], K: int = 0) -> list:
    """
    Load Q files for alignment within K clusters.
    
    Args:
        qfiles: List of Q file paths.
        K: Number of clusters.
        
    Returns:
        List of matrices loaded from the Q files.
    """
    if not qfiles:
        raise ValueError("No Q files provided for loading.")
    
    if K == 0:
        logger.warning("K was not specified, using the number of columns from the first file.")

    Q_list = []
    for file_path in qfiles:
        try:
            matrix = load_matrix(file_path)
            if K == 0:
                K = matrix.shape[1]
            if matrix.shape[1] != K:
                raise ValueError(f"File {file_path} does not have {K} columns (clusters).")
            Q_list.append(matrix)
        except Exception as e:
            logger.error(f"ERROR loading {file_path}: {e}")
    
    if K == 0:
        logger.error("No successful file reading.")
    else:
        logger.info(f"Successfully loaded {len(Q_list)} Q files with K={K} clusters.")

    return Q_list


def write_aligned_within_k(aligned_summary: list, Q_names: list, output_file: str):
    """
    Write within-K aligned Q results to output files.
    
    Args:
        aligned_summary: List of aligned result summaries.
        Q_names: List of Q file base names (without extensions).
        output_file: Directory to save the output result file.
    """
    
    with open(output_file,"w") as f:
        f.write('Replicate1-Replicate2,Cost,Alignment\n')
        for i_row in range(len(aligned_summary)):
            i,j = aligned_summary[i_row][0]
            if not (isinstance(i, int) and isinstance(j, int)):
                raise ValueError("The first element of each item in aligned_summary must be a tuple of two integers (i, j).")
            id_i = Q_names[i]
            id_j = Q_names[j]
            cost = float(aligned_summary[i_row][1])
            pattern = pattern_to_str(aligned_summary[i_row][2])
            pattern_rev = pattern_to_str(reverse_alignment_same_K(aligned_summary[i_row][2]))
            f.write(f"{id_i}-{id_j},{cost},{pattern}\n")
            f.write(f"{id_j}-{id_i},{cost},{pattern_rev}\n")


def load_aligned_within_k(input_file: str, Q_names: List[str] = list()) -> list:
    """
    Load within-K aligned Q results from an input file.
    
    Args:
        input_file: Path to the input result file.
        Q_names: Optional list of Q file base names (without extensions) to map names to indices.
        
    Returns:
        List of aligned result summaries.
    """
    aligned_summary = []
    with open(input_file, "r") as f:
        header = f.readline().strip().split(',')
        if header != ['Replicate1-Replicate2', 'Cost', 'Alignment']:
            raise ValueError("Input file format is incorrect or header is missing.")
        
        for line in f:
            parts = line.strip().split(',')
            if len(parts) != 3:
                logger.warning(f"Skipping malformed line: {line.strip()}")
                continue
            replicate_pair, cost_str, pattern_str = parts
            try:
                id1, id2 = replicate_pair.split('-')
                if not Q_names:
                    i, j = id1, id2
                else:
                    i = Q_names.index(id1)
                    j = Q_names.index(id2)
                cost = float(cost_str)
                pattern = str_to_pattern(pattern_str)
                aligned_summary.append([(i, j), cost, pattern])
            except Exception as e:
                logger.error(f"Error processing line '{line.strip()}': {e}")
    
    return aligned_summary


def aligned_res_to_dicts(aligned_summary: List[Any], Q_names: List[str]) -> Tuple[dict, dict]:
    """ Convert aligned summary list to a dictionary for easy lookup.
        Args:
            aligned_summary: List of aligned result summaries.
            Q_names: List of Q file base names (without extensions).
        Returns:
            Dictionary mapping pairs of Q names to their alignment patterns (both ordering).
            Dictionary mapping pairs of Q names to their alignment costs (i<j).
    """
    alignment_withinK = dict()
    cost_withinK = dict()
    for i_row in range(len(aligned_summary)):
        i,j = aligned_summary[i_row][0]
        id_i = Q_names[i]
        id_j = Q_names[j]
        cost = float(aligned_summary[i_row][1])
        pattern = aligned_summary[i_row][2]
        alignment_withinK[(id_i,id_j)] = pattern
        alignment_withinK[(id_j,id_i)] = reverse_alignment_same_K(pattern)
        cost_withinK[(id_i,id_j)] = cost
    return alignment_withinK, cost_withinK