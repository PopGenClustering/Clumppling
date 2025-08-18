import numpy as np
from typing import List

import logging  
logger = logging.getLogger(__name__)

def cd_custom(adj_mat: np.ndarray) -> List[int]:
    
    n_nodes = adj_mat.shape[0]

    # Please comment out the following line and customize your community detection method here.
    # ## Dummuy placeholder
    # communities = list(range(n_nodes))  
    # logger.warning(f"Using dummy custom community detection method. Output will be each node in its own community (mode).")

    # ## Example 1: use Markov clustering
    import markov_clustering
    # Run MCL
    result = markov_clustering.run_mcl(adj_mat, pruning_threshold=0.2)
    coms = markov_clustering.get_clusters(result)
    communities = [-1] * n_nodes
    for cluster_idx, nodes in enumerate(coms):
        for node in nodes:
            communities[node] = cluster_idx

    ## Example 2: use infomap
    # import networkx as nx
    # from cdlib import algorithms
    # G = nx.from_numpy_array(adj_mat*100)  # Scale the adjacency matrix
    # G.remove_edges_from(nx.selfloop_edges(G))
    # coms = algorithms.infomap(G)
    # logger.info(f"Custom method (example 2. infomap clustering with adjacency matrix scaled by 100) detected {len(coms.communities)} communities.")
    # communities = [-1] * n_nodes
    # for cluster_idx, nodes in enumerate(coms.communities):
    #     for node in nodes:
    #         communities[node] = cluster_idx

    return communities
