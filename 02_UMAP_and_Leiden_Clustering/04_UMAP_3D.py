# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 20:40:07 2021
@author: ezunder
"""

import numpy as np
import pandas as pd
import umap
import hnswlib

# Input/output file names
INPUT_EXPRS = "expression_matrix.csv"
INPUT_PANEL = "panel.csv"
UMAP_XYZ_OUTPUT = "umap_xyz_3D.csv"

# Parameters for HNSW approximate nearest neighbor determination
# Currently set to default values.  Not sure how changing these will effect results or runtime. . .
HNSW_SPACE = 'l2' #'l2' Squared L2: (sum((Ai-Bi)^2)), 'ip' Inner Product: (1.0 - sum(Ai*Bi)), or 'cosine' Cosine Similarity: 1.0 - sum(Ai*Bi) / sqrt(sum(Ai*Ai) * sum(Bi*Bi))
HNSW_EF_CONSTRUCTION = 200
HNSW_M = 16
HNSW_SET_EF = 20

# Parameters for UMAP embedding
UMAP_N_NEIGHBORS = 15 # lower=local emphasis, higher=global emphasis
UMAP_RANDOM_SEED = 99 # See how much it changes from run to run!
UMAP_MIN_DIST = 0
UMAP_N_COMPONENTS = 3 # Can try 3D if 2D has too much overlap
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

# Convert to numpy array for input to UMAP
umap_np = np.array(pd.read_csv(INPUT_EXPRS))

def run_umap(umap_np):

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
    )

    ## Get the optimal a and b embedding values
    ## (regular UMAP function umap.umap_.UMAP() takes care of this automatically,
    ## but simplicial_set_embedding doesn't, so we do it here instead)
    a, b = umap.umap_.find_ab_params(UMAP_SPREAD, UMAP_MIN_DIST)

    ## Get the simplicial set embedding (this is like running UMAP on our own graph)
    data_embed = umap.umap_.simplicial_set_embedding(
        data=umap_np,
        graph=fuzzy_set[0],
        n_components=UMAP_N_COMPONENTS,
        initial_alpha=UMAP_INITIAL_ALPHA,
        a=a,
        b=b,
        gamma=UMAP_GAMMA,
        negative_sample_rate=UMAP_NEGATIVE_SAMPLE_RATE,
        n_epochs=UMAP_N_EPOCHS,
        init=UMAP_INIT,
        random_state=np.random.RandomState(seed=UMAP_RANDOM_SEED),
        metric=UMAP_METRIC,
        metric_kwds=UMAP_METRIC_KWDS,
        densmap=UMAP_DENSMAP,
        densmap_kwds=UMAP_DENSMAP_KWDS,
        output_dens=UMAP_OUTPUT_DENS,
        output_metric=UMAP_OUTPUT_METRIC,
        output_metric_kwds=UMAP_OUTPUT_METRIC_KWDS,
        euclidean_output=UMAP_EUCLIDEAN_OUTPUT,
        parallel=UMAP_PARALLEL,
        verbose=UMAP_VERBOSE)[0] # Outputs more than just sparse matrix, so just take first element (sparse matrix)

    # ## This is the standard way of running umap, but it takes longer because it
    # ## doesn't let you import your own knn values.
    # data_embed = umap.UMAP(n_neighbors=UMAP_N_NEIGHBORS, min_dist=UMAP_MIN_DIST,
    #                  n_components=UMAP_N_COMPONENTS, n_epochs=UMAP_N_EPOCHS,
    #                  metric=UMAP_METRIC, ).fit_transform(umap_np)
    return data_embed
  
data_embed = run_umap(umap_np)
# Extract x and y coordinates from UMAP output
out_df = pd.DataFrame({'umap_x': data_embed[:, 0],'umap_y': data_embed[:, 1], 'umap_z': data_embed[:,2]})

# Output umap XYZ coordinates
out_df.to_csv(UMAP_XYZ_OUTPUT, index=False)


