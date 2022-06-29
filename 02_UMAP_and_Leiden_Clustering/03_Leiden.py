# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 04:25:48 2021
@author: Eli Zunder
"""
## Import packages
import pandas as pd
import numpy as np
import leidenalg
import igraph as ig
#import csv
#import os
import umap
import hnswlib
#import time

# Need to change this for the level of subsetting/subclustering!
CLUSTER_ROUND = 1 # start with 1, then do 2 for subclustering, then do 3 for sub-subclustering, etc.

# Input/output file names
CLUSTER_BASENAME = "cluster_R" # For output files
ASSIGNS_BASENAME = "_assigns.csv" # For cluster assignment file
GROUPS_BASENAME = "_groups.csv" # For cluster grouping file
INPUT_EXPRS = "expression_matrix_analysis.csv" # Expression matrix dataset
INPUT_PANEL = "panel.csv" # Panel info, for which markers to use, etc.

# Parameters for HNSW approximate nearest neighbor determination
# Currently set to default values.  Not sure how changing these will effect results or runtime. . .
HNSW_SPACE = 'l2' #'l2' Squared L2: (sum((Ai-Bi)^2)), 'ip' Inner Product: (1.0 - sum(Ai*Bi)), or 'cosine' Cosine Similarity: 1.0 - sum(Ai*Bi) / sqrt(sum(Ai*Ai) * sum(Bi*Bi))
HNSW_EF_CONSTRUCTION = 200
HNSW_M = 16
HNSW_SET_EF = 20

# Parameters for UMAP embedding
UMAP_N_NEIGHBORS = 15 # lower=local emphasis, higher=global emphasis
UMAP_RANDOM_SEED = 42 # See how much it changes from run to run!
UMAP_MIN_DIST = 0
UMAP_N_COMPONENTS = 2 # Can try 3D if 2D has too much overlap
UMAP_N_EPOCHS = 1000 # Does running longer give a better-defined layout?
UMAP_METRIC = 'euclidean' # Cosine distance is less scale-dependent
UMAP_ANGULAR = False
UMAP_SET_OP_MIX_RATIO = 1.0
UMAP_LOCAL_CONNECTIVITY = 1.0
UMAP_VERBOSE = True # Might want to turn this off some times?
UMAP_SPREAD = 1.0 # Lowest spread (~1.0) generally looks "better" and less blobby
UMAP_INITIAL_ALPHA = 1.0
UMAP_GAMMA = 1.0
UMAP_NEGATIVE_SAMPLE_RATE = 5
UMAP_N_EPOCHS = 1000
UMAP_INIT = 'spectral'
UMAP_METRIC_KWDS = {}
UMAP_DENSMAP = False
UMAP_DENSMAP_KWDS = {}
UMAP_OUTPUT_DENS = False
UMAP_OUTPUT_METRIC = 'euclidean'
UMAP_OUTPUT_METRIC_KWDS = {}
UMAP_EUCLIDEAN_OUTPUT = True
UMAP_PARALLEL = False

UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0 # Like prevent-overlap in Gephi
UMAP_N_COMPONENTS = 2 # Dimensions for embedding/layout (usually want 2D or 3D)
UMAP_N_EPOCHS = 1000 # How long to let the force-directed layout/gradient descent to run
UMAP_METRIC = 'euclidean' # Many other metrics possible.  Manhattan, Cosine, etc.

def fuzzy_leiden(umap_np):

    ## Use HNSW to find approximate nearest neighbors (ANN)
    p = hnswlib.Index(space = HNSW_SPACE, dim = umap_np.shape[1])
    data_labels = np.arange(umap_np.shape[0])
    p.init_index(max_elements = umap_np.shape[0],
                 ef_construction = HNSW_EF_CONSTRUCTION, M = HNSW_M)
    p.add_items(umap_np, data_labels)
    p.set_ef(HNSW_SET_EF) # ef should always be > k
    labels, distances = p.knn_query(umap_np, k = UMAP_N_NEIGHBORS)

    ## Convert ANN to fuzzy simplicial set (think of this like the knn UMAP uses)
    fuzzy_set = umap.umap_.fuzzy_simplicial_set(
        X=umap_np,
        n_neighbors=UMAP_N_NEIGHBORS,
        random_state=np.random.RandomState(seed=UMAP_RANDOM_SEED),
        metric=UMAP_METRIC,
        metric_kwds={},
        knn_indices=labels,
        knn_dists=distances,
        angular=UMAP_ANGULAR,
        set_op_mix_ratio=UMAP_SET_OP_MIX_RATIO,
        local_connectivity=UMAP_LOCAL_CONNECTIVITY,
        verbose=UMAP_VERBOSE,
    )[0] # Just take the sparse matrix, not the extras

    ## Build igraph graph from fuzzy simplicial set sparse matrix
    leiden_graph = ig.Graph(np.column_stack((fuzzy_set.nonzero())), directed=False)
    weights = np.array(1-fuzzy_set.data) # do 1-x instead of 1/x because range is from 0-1
    leiden_graph.es['weight'] = weights

    ## Do Leiden clustering on fuzzy simplicial set graph
    leiden_out = leidenalg.find_partition(leiden_graph,
                                          leidenalg.ModularityVertexPartition,
                                          n_iterations=-1,seed=1)

    # Add one to every cluster ID to make it 1-indexed instead of 0-indexed!
    plus_one = [x+1 for x in leiden_out.membership]
    return plus_one

# Convert to numpy array for input to UMAP
umap_np = np.array(pd.read_csv(INPUT_EXPRS))

if CLUSTER_ROUND > 1:
    # If round > 1, then use previous round clustering from cluster_ids file
    # and grouping from cluster metadata file to subset, and then do Leiden
    # clustering on each subset.
    assigns_infile = CLUSTER_BASENAME + str(CLUSTER_ROUND - 1) + ASSIGNS_BASENAME
    groups_infile = CLUSTER_BASENAME + str(CLUSTER_ROUND - 1) + GROUPS_BASENAME

    assigns_in = pd.read_csv(assigns_infile)
    assigns_prev = assigns_in['R' + str(CLUSTER_ROUND - 1)]
    groups_in = pd.read_csv(groups_infile)

    out_clusts = [0]*len(umap_np)
    add_for_unique = 0 #use this to clusters from separate subclustering rounds unique ids

    # Carry through 0-labeled clusters (don't subcluster them).
    # Rename them [0:n] to keep everything in order nicely.
    zero_clusts = np.array(groups_in.loc[groups_in['Group_Assigns'] == 0, 'R' + str(CLUSTER_ROUND - 1)])
    if len(zero_clusts) > 0:
        for c in range(len(zero_clusts)):
            for i in np.where(assigns_prev==zero_clusts[c])[0]:
                out_clusts[i] = c+1
        add_for_unique = c+1

    # Get all non-zero groups.  (Already took care of the zero group by carrying it through)
    all_groups = np.unique(groups_in['Group_Assigns'])
    nonzero_groups = all_groups[all_groups > 0]

    for g in nonzero_groups:
        g_clusts = np.array(groups_in.loc[groups_in['Group_Assigns'] == g, 'R' + str(CLUSTER_ROUND - 1)])
        umap_np_sub = umap_np[assigns_in['R' + str(CLUSTER_ROUND - 1)].isin(g_clusts)]
        clust_ids = fuzzy_leiden(umap_np_sub)
        clust_ids = [x+add_for_unique for x in clust_ids]
        add_for_unique = np.max(clust_ids)
        replace_idxs = np.where(assigns_in['R' + str(CLUSTER_ROUND - 1)].isin(g_clusts))[0].tolist()
        for (replace_idxs, clust_ids) in zip(replace_idxs, clust_ids):
            out_clusts[replace_idxs] = clust_ids
else:
    # If round 1, then just do regular clustering
    out_clusts = fuzzy_leiden(umap_np)


## Output cluster ids
assigns_outfile = CLUSTER_BASENAME + str(CLUSTER_ROUND) + ASSIGNS_BASENAME
assigns_out = pd.DataFrame(out_clusts, columns=['R' + str(CLUSTER_ROUND)])
assigns_out.to_csv(assigns_outfile, encoding='utf-8', index=False)

## Output cluster metadata - this file will have to be manually updated for cluster grouping
groups_outfile = CLUSTER_BASENAME + str(CLUSTER_ROUND) + GROUPS_BASENAME
group_vals = np.unique(out_clusts)
groups_out = pd.DataFrame(group_vals, columns=['R' + str(CLUSTER_ROUND)])
groups_fill_cols = pd.DataFrame(np.full([len(group_vals),4], ""),
                               columns=['Cluster_Labels', 'Group_Assigns',
                                        'Groups_Unique', 'Group_Labels'])
groups_out = pd.concat([groups_out, groups_fill_cols], axis=1)
#groups_out = groups_out.replace(np.nan, '', regex=True)
groups_out.to_csv(groups_outfile, encoding='utf-8', index=False)
