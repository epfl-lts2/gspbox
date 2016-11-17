Matlab / C++ implementation of community detection algorithm.
After "Fast unfolding of community hierarchies in large networks"
Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte and
Etienne Lefebvre
Journal of Statistical Mechanics: Theory and Experiment, 1742-5468, P10008 (12 pp.), 2008.

Implementation : Antoine Scherrer
antoine.scherrer@ens-lyon.fr

*** USAGE - Full matlab

See help of m files for details.

cluster_jl.m : Weighted (or not), non oriented version of algorithm 
 matrix is symetrized using sum of incoming and outgoing weights)

cluster_jl_orient.m : Weighted (or not), oriented version of algorithm 
 using extended definition of modularity for oriented graphs 

cluster_jl_orientT.m : Weighted (or not), oriented version of algorithm 
 using symetric matrix A = M*M^t 

*** USAGE - Matlab/C++

You need to compile jl_clust.cpp, jl_mnew.cpp and jl_clust_orient.cpp
with mex compiler. Then you can use the following routines to perform
the community detection faster.

cluster_jl_cpp.m : Weighted (or not), non oriented version of algorithm 
 matrix is symetrized using sum of incoming and outgoing weights)

cluster_jl_orient_cpp.m : Weighted (or not), oriented version of algorithm 
 using extended definition of modularity for oriented graphs 

cluster_jl_orientT_cpp.m : Weighted (or not), oriented version of algorithm 
 using symetric matrix A = M*M^t 


