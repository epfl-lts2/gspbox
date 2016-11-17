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
%   gsp_sort_nodes(W) sorts the nodes of the weighted adjacency matrix W
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
%   Example:
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
%   References:
%     V. D. Blondel, J.-L. Guillaume, R. Lambiotte, and E. Lefebvre. Fast
%     unfolding of communities in large networks. Journal of statistical
%     mechanics: theory and experiment, 2008(10):P10008, 2008.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_sort_nodes.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

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



