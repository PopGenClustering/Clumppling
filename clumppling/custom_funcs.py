# -*- coding: utf-8 -*-
"""
Clumppling: customizable functions

@author: Xiran Liu
"""

import networkx as nx
import community.community_louvain as community_louvain

def cd_custom(G):
    """Customized community detection method (need to be modified)

    Parameters
    ----------
    G : networkx.Graph
        the similarity network of replicates, with edges weighted by similarity after optimal alignment

    Returns
    -------
    partition_map
        a dictionary where keys are the indices of replicates and values are the indices of the communities they belong to
    
    Please comment out the dummy placeholder and customize your community detection method there.
    Two examples are provided below in commented text.
    """

    # Please comment out the following line and customize your community detection method here.
    ## Dummuy placeholder
    partition_map = {i:0 for i in range(G.number_of_nodes())} 

    ## Example 1: use Markov clustering
    # from cdlib import algorithms
    # coms = algorithms.markov_clustering(G,pruning_threshold=0.2) # scipy 1.8.0
    # partition = coms.communities
    # partition_map = {i:i_s for i_s,s in enumerate(partition) for i in s}

    ## Example 2: use Louvain with additional modularity check
    # resolution = 1.0 
    # partition_map = community_louvain.best_partition(G,resolution=resolution,random_state=6)
    # cd_mod_thre = 0.2
    # mod_res = community_louvain.modularity(partition_map, G)
    # if mod_res<cd_mod_thre: # quality of community detection is low --> no community structure
    #     partition_map = {i:0 for i in range(G.number_of_nodes())} 
    
    return partition_map