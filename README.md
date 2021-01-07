# CIRCA Clustering
Here we present a clustering function that finds the optimal clustering of an input similarity matrix by maximizing the Rand index (RI) via coordinate ascent.
In short, a proposal is made to move a random stimulus (i.e., row/column) to a random, different cluster.
If this proposal generates an increase in the RI, the proposal is accepted.
This procedure is executed iteratively until either the desired number of coordinate ascent iterations is met, or until the algorithm converges, and no improvements can be made.

CIRCA_Clustering.m is a MATLAB function that accepts three arguments:
1. hedges - a symmetric n x n (where n is the number of stimuli) similarity/adjacency matrix that codes the pairwise similarity between stimuli in the design. Missing values must be represented as NaN values, and all values must be between 0 and 1. 
2. ng - number of groups/clusters/categories
3. maxits - number of coordinate ascent iterations. To ensure convergence, set this to inf. 

And there are 5 outputs:
1. gs - An nx1 vector of integers from 1-ng, which gives the cluster assignments for every row/column (i.e., image/item) in hedges.
2. ngcount - number of moves made via coordinate ascent.
3. counter - number of proposals made. 
4. R - Rand Index between similarity matrix and optimized clustering.
5. converged - true if the algorithm reaches a stationary point (and no moves that improve the rand index are possible)

Example.m is a MATLAB script that shows how CIRCA_Clustering.m works on a simulated dataset, with a small amount of response noise. 

See our paper (Category Systems for Real World Scenes, Anderson et al. (2021)) for further details. 
