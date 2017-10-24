function [W_sorted, ind_sort, modularity, com_id, com_sizes] = gsp_sort_nodes(W, recursive, self_loops)
%GSP_SORT_NODES Sort nodes using the Louvain clustering method
%   Usage:  [Wout, ind, modularity, com_id, com_sizes] = gsp_sort_nodes(W)
%           [Wout, ind, modularity, com_id, com_sizes] = gsp_sort_nodes(W, recursive)
%           [Wout, ind, modularity, com_id, com_sizes] = gsp_sort_nodes(W, recursive, self_loops)
%
%   Inputs:
%         W         : Weigthed adjacency matrix. Can contain self-loops
%         recursive : Use recursive computation (default: 1)
%         self_loops: Take into account self loops (default: 0)
%
%   Outputs:
%         Wout      : Sorted weighted adjacency matrix
%         ind       : Indices used for sorting:       Wout = W(ind, ind)
%         modularity: Modularity of final clustering used for sorting
%         com_id    : Community ID of each node in the initial W
%         com_sizes : Number of nodes in each community
%
%   gsp_sort_nodes(W) sorts the nodes of the weighted adjacency matrix $W$
%   by clustering them according to the "Louvain" method.
%
%   By default the clustering is recursive and the level of recursion that
%   gives the maximum modularity is kept. The nodes of the graph are sorted
%   according to the clusters given.
%
%   The clustering method is using code given free online by Antoine
%   Scherrer available here: 
%   https://perso.uclouvain.be/vincent.blondel/research/louvain.html 
%
%   Example:::
%
%         G = gsp_nn_graph([[randn(100, 1)-1; (randn(100, 1)+7)], ...
%           [randn(50, 1); (randn(100,1)-4); (randn(50, 1)*2+8)]]);
%         [Wout, ind, modularity, com_id, com_sizes] = gsp_sort_nodes(G.W);
%         figure; imagesc(Wout); title('sorted adjacency matrix');
%         figure; gsp_plot_signal(G, ind);
%       
%
%
%   See also: cluster_jl
% 
%   References: blondel2008fast
%

% Author: Vassilis Kalofolias
% Date: June 2016


if nargin < 2
    % Use recursive clustering and pick the level with the highest modularity
    recursive = 1;
end
if nargin < 3
    % Don't use self loops by default
    self_loops = 0;
end
debug = 0;
verbose = 0;



COMTY = cluster_jl(W, recursive, self_loops, debug, verbose);

% Keep as best the splitting with highest modularity
ind_best_splitting = argmax(COMTY.MOD);

% Which node belongs to which community
com_id = COMTY.COM{ind_best_splitting};

% Sort according to the communities
[~, ind_sort] = sort(com_id);

modularity = COMTY.MOD(ind_best_splitting);
com_sizes = COMTY.SIZE{ind_best_splitting};

W_sorted = W(ind_sort, ind_sort);


