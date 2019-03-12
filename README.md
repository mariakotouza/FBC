# FBC
A multi-metric frequency based hierarchical clustering framework
The framework consists of 3 major modules:
1) The binary tree construction algorithm that creates a hierarchy of the documents using three metrics (Identity, Entropy, Bin Similarity).
There are two available approaches to create construct the binary tree: per cluster (binary_tree_construction_per_cluster.R) and per level (binary_tree_construction_per_level.R).4
The first one is prefered for smaller datasets, whereas the second one is quicer for larger dataasets and can run in parallel using different cores.
2) The branch breaking algorithm which composes the final clusters by applying thresholds to each branch of the tree (branch_breaking.R).
3) The clustering algorithm is followed by a meta-clustering module which makes use of the graph theory to obtain insights in the leaf clusters' connections. (distances_and_graph.R)