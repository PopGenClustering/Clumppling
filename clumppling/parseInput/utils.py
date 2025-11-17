import os
import numpy as np
from collections import OrderedDict
from typing import List, Optional, Tuple
import logging
logger = logging.getLogger(__name__)

def group_labels(labels: list, ordered_uniq_labels: Optional[List[str]] = None) -> Tuple[list, list]:
    """ Group identical labels together while preserving the order they first appear.
        Args:
            labels: List of labels.
        Returns:
            grouped_labels: List of labels with identical labels grouped together.
            reorder_indices: List of indices to reorder the original labels.
    """
    # Group labels while preserving first-seen order, and if ordered_uniq_labels provided, use that order
    if ordered_uniq_labels:
        # check if ordered_uniq_labels match labels
        assert set(ordered_uniq_labels) == set(labels), "ordered_uniq_labels does not match labels."
        assert len(ordered_uniq_labels) == len(set(labels)), "ordered_uniq_labels contains duplicates."
        groups = OrderedDict((label, []) for label in ordered_uniq_labels)
    else:
        groups = OrderedDict()
    for idx, label in enumerate(labels):
        groups.setdefault(label, []).append(idx)  # store indices instead of label

    # Flatten indices to get the reordering
    reorder_indices = [i for idx_list in groups.values() for i in idx_list]
    # Apply to labels to get grouped labels
    grouped_labels = [labels[i] for i in reorder_indices]
    logger.info(f"Out of {len(labels)} labels, {len(groups)} are unique.")

    return grouped_labels, reorder_indices


def process_files(
    input_dir: str, output_dir: str,
    fmt: str = "generalQ", extension: str = "", 
    skip_missing: bool = True, 
    delimiter: str = " ", skip_rows: int = 0,
    tolerance: float = 1e-6, meta_file: str = "input_meta.txt",
    label_cols: list[int]=[0, 1, 3], mat_start_col: int=5,
    reorder_indices: Optional[List[int]] = None,
    reorder_uniq_labels: Optional[List[str]] = None
) -> Optional[np.ndarray]:
    
    file_paths = [os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith(extension)]
    if not file_paths:
        logger.error(f"No files found with the specified extension '{extension}'.")
        return
    
    grouped_files = {}
    metadata = []
    expected_nrows = None
    labels = None

    for file_path in file_paths:
        try:
            if fmt =="admixture" and extension == ".indivq":
                labels, matrix = extract_labels_and_matrix(
                    file_path, skip_missing, delimiter, skip_rows,
                    label_cols=[0, 1, 3], mat_start_col=5
                )
                if labels is not None:
                    labels = [label[-1] for label in labels]
                if labels is None or matrix is None:
                    logger.error(f"Failed to extract labels or matrix from {file_path}.")
                    continue
            elif fmt == "structure":
                tolerance = 1e-2  # Default tolerance for STRUCTURE files
                if extension != "_f":
                    logger.warning(f"STRUCTURE file extension is not '_f' for {file_path}.")
                lines = extract_lines_from_structure_f(file_path)
                labels, matrix = extract_labels_and_matrix_from_lines(
                    lines, skip_missing, delimiter, skip_rows,
                    label_cols=label_cols, mat_start_col=mat_start_col
                )
                if labels is not None:
                    labels = [label[-1] for label in labels]  # Use the last column as label
                if labels is None or matrix is None:
                    logger.error(f"Failed to extract labels or matrix from {file_path}.")
                    continue
            elif fmt in ["generalQ", "admixture", "fastStructure"]:
                matrix = load_matrix(file_path, skip_missing, delimiter, skip_rows)
            else:
                logger.error(f"Unsupported format '{fmt}' for file {file_path}.")
                continue
            
            validate_membership_matrix(matrix, tolerance)
            matrix = rescale_membership_matrix(matrix)

            # Check number of columns (clusters)
            n_rows, n_cols = matrix.shape
            if n_cols == 1:
                logger.warning(f"Ignored file {file_path} with K=1.")
                continue

            # Check number of rows (individuals) consistency
            if expected_nrows is None:
                expected_nrows = n_rows
            elif n_rows != expected_nrows:
                logger.error(f"Row number mismatch in {file_path} (expected {expected_nrows}, got {n_rows}). Skipping.")
                continue
            
            # Group and store
            grouped_files.setdefault(n_cols, []).append((file_path, matrix))

        except ValueError as e:
            logger.error(f"ERROR processing {file_path}: {e}")
    
    # Save labels if applicable
    if labels is not None:        
        if reorder_indices is None:
            # reorder    
            labels, reorder_indices = group_labels(list(labels), ordered_uniq_labels=reorder_uniq_labels)
            labels = np.array(labels)
            if not reorder_indices == list(range(len(labels))):
                logger.warning(f"Individual labels reordered by grouping identical labels together.")

        labels_out_path = os.path.join(output_dir, "input_labels.txt")
        np.savetxt(labels_out_path, labels, fmt="%s", delimiter=',')
        logger.info(f"Labels saved to {labels_out_path}")
    

    global_id = 1  # ID
    for k in sorted(grouped_files.keys()):  # Loop in increasing K
        group = grouped_files[k]
        logger.info(f"Found {len(group)} files with K={k}")

        for i, (file_path, matrix) in enumerate(group):
            filename = f"{global_id}_K{k}R{i+1}.Q"
            out_path = os.path.join(output_dir, filename)
            if reorder_indices is not None:
                matrix = matrix[reorder_indices, :]
            np.savetxt(out_path, matrix, delimiter=' ') #fmt="%.6f",
            metadata.append(f"{file_path},{out_path},{k}")
            global_id += 1

    # Save metadata
    meta_out_path = os.path.join(output_dir,meta_file) 
    with open(meta_out_path, "w") as f:
        f.write("\n".join(metadata))
    logger.info(f"Metadata written to {meta_out_path}")

    return labels

def load_matrix(file_path: str, skip_missing: bool = True, 
                delimiter: str = " ", skip_rows: int = 0):
    matrix = []
    expected_ncols = None
    if skip_rows > 0:
        logger.info(f"Skipping the first {skip_rows} rows in {file_path}.")
    
    with open(file_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if line_num < skip_rows:
                continue
            line = line.strip()
            if not line:
                continue
            entries = line.split() if delimiter == " " else line.split(delimiter)
            if expected_ncols is None:
                expected_ncols = len(entries)
            elif (len(entries) != expected_ncols):
                if skip_missing:
                    logger.warning(f"Skipping line {line_num} in {file_path} due to missing data.")
                    continue
                else:
                    raise ValueError(f"Missing data at line {line_num}.")
            try:
                row = list(map(float, entries))
                matrix.append(row)
            except ValueError:
                raise ValueError(f"Non-numeric entry at line {line_num}.")
    if not matrix:
        raise ValueError("File is empty or all lines were skipped.")
    return np.array(matrix)


def extract_labels_and_matrix_from_lines(lines: List[str], skip_missing: bool = True,
                                         delimiter: str = " ", skip_rows: int = 0,
                                         label_cols: Optional[List[int]] = None,
                                         mat_start_col: int = 2) -> Tuple[Optional[np.ndarray], np.ndarray]:
    labels = []
    matrix_rows = []
    expected_ncols = None

    for line_num, line in enumerate(lines):
        if line_num < skip_rows:
            continue

        line = line.strip()
        if not line:
            continue  # skip empty lines

        # Example line tokens:
        # ['1', 'HGDP00904', '(0)', '1', ':', '0.000010', '0.999990']
        parts = line.split() if delimiter == " " else line.split(delimiter)
        try:
            # Extract float values in matrix part
            entries = list(map(float, parts[mat_start_col:]))
            if expected_ncols is None:
                expected_ncols = len(entries)
            elif (len(entries) != expected_ncols):
                if skip_missing:
                    logger.error(f"Column number mismatch in line {line_num} (expected {expected_ncols}, got {len(entries)}). Skipping.")
                    continue
                else:
                    raise ValueError(f"Missing data at line {line_num}.")
            # Extract labels if requested
            if label_cols is not None:
                label = [parts[idx_col] for idx_col in label_cols]
                labels.append(label)
            matrix_rows.append(entries)

        except Exception as e:
            logger.warning(f"Error parsing line {line_num}: {e}")

    if not matrix_rows:
        raise ValueError("No valid data found in file.")
    
    labels = np.array(labels) if len(labels)>0 else None
    matrix = np.array(matrix_rows)

    return labels, matrix


def extract_labels_and_matrix(file_path: str, skip_missing: bool = True,
                              delimiter: str = " ", skip_rows: int = 0,
                              label_cols: Optional[List[int]] = None, 
                              mat_start_col: int = 2) -> Tuple[Optional[np.ndarray], np.ndarray]:
    
    if label_cols is None and mat_start_col == 0:
        return None, load_matrix(file_path, skip_missing, delimiter, skip_rows)

    with open(file_path, "r") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"File {file_path} is empty or all lines were skipped.")

    labels, matrix = extract_labels_and_matrix_from_lines(
        lines, skip_missing, delimiter, skip_rows,
        label_cols=label_cols, mat_start_col=mat_start_col
    )

    return labels, matrix
    

def extract_lines_from_structure_f(file_path: str) -> List[str]:
    """
    Extracts lines with labels and membership matrix from a STRUCTURE file (_f).
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"File {file_path} is empty or all lines were skipped.")
            
    # get only the membership coefficients part
    res_start = lines.index("Inferred ancestry of individuals:\n")
    lines = lines[(res_start+2):]
    res_end = lines.index("\n")
    lines = lines[:res_end]
            
    return lines


def validate_membership_matrix(matrix: np.ndarray, tolerance: float =  1e-6):
    """ 
    Validates the membership matrix to ensure all values are in [0, 1] 
    and each row sums to 1 within a tolerance. 
    """
    if not ((0 <= matrix) & (matrix <= 1)).all():
        raise ValueError("Matrix contains values outside [0, 1].")

    row_sums = matrix.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=tolerance):
        bad_rows = np.where(~np.isclose(row_sums, 1.0, atol=tolerance))[0]
        raise ValueError(f"Rows {bad_rows.tolist()} do not sum to 1 (±{tolerance}).")
    

def rescale_membership_matrix(matrix: np.ndarray, tolerance: float = 1e-9) -> np.ndarray:
    """
    Rescales the membership matrix to ensure each row sums to 1.
    """
    # add a small value to avoid division by zero
    matrix += tolerance
    row_sums = matrix.sum(axis=1)

    return matrix / row_sums[:, np.newaxis]


def extract_meta_input(output_dir: str, meta_file: str = "input_meta.txt") -> Tuple[List[str], List[int], dict]:
    """ Extracts metadata from the input files processed by `process_files`.
    Args:
        output_dir: Directory where the processed files are saved.
        meta_file: Name of the metadata file (default: input_meta.txt).
    Returns:
        A tuple containing:
        - A list of Q file names (without paths).
        - A list of K values corresponding to each Q file. 
        - A dictionary mapping each K to a list of individual IDs.
    """     
    meta_out_path = os.path.join(output_dir,meta_file) 
    if not os.path.exists(meta_out_path):
        raise FileNotFoundError(f"Metadata file {meta_out_path} does not exist.")
    meta_data = np.loadtxt(meta_out_path, dtype=str, delimiter=',')
    K_list = meta_data[:,2].astype(int).tolist()
    K_range = np.sort(np.unique(K_list)).tolist()
    Q_files = meta_data[:,1].tolist()
    Q_names = [os.path.splitext(os.path.basename(p))[0] for p in Q_files]
    IDs = [int(p.split('_')[0]) for p in Q_names]  # Extract IDs from file names
    K2IDs = {k:[] for k in K_range}
    for i_r,r in enumerate(IDs):
        K = K_list[i_r]
        K2IDs[K].append(r)

    return Q_names, K_range, K2IDs

    