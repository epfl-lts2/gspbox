function [G] = gsp_adj2vec(G)
%GSP_ADJ2VEC Prepare the graph for the gradient computation
%   Usage: [G] = gsp_adj2vec(G)
%
%   Input parameters:
%       G   : Graph structure
%   Output parameters:
%       G   : Graph structure
%
%   This function converts adjacency matrix to edge vector form. It also
%   add the field G.Diff that is the sparse gradient matrix
%
%   See also: gsp_grad gsp_div
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_adj2vec.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

% Author: Nathanael Perraudin, Vassilis Kalofolias
% Date  : 5 August 2014
% Testing: test_operators

if isfield(G,'Gm');
    G.Gm = gsp_adj2vec(G.Gm);
end

if G.directed
    error('GSP_ADJ2VEC: Not implemented yet');
else
    % Keep each edge only once (they are duplicated!). Keep also loops.
    if G.directed
        error('Not implemented now!')
    end
    [v_i, v_j, weights] = find(tril(G.W));
    G.v_in = v_i;
    G.v_out = v_j;

    % the indices of the edges in the Adgacency matrix A:
    %G.ind_edges = sub2ind(size(G.W), G.v_in, G.v_out);
    G.weights = weights;         % |E| x 1
    G.Ne = length(G.v_in);

    G.Diff = gsp_grad_mat( G );
end

end




% Old code from Vasilis
% [v_i, v_j] = find(G.W);
% 
% % Keep each edge only once (they are duplicated!). Keep also loops.
% ind_keep = (v_i >= v_j);
% G.v_in = v_i(ind_keep);
% G.v_out = v_j(ind_keep);
% 
% % the indices of the edges in the Adgacency matrix A:
% G.ind_edges = sub2ind(size(G.W), G.v_in, G.v_out);
% G.weights = full(G.W(G.ind_edges));         % |E| x 1
% G.Ne = size(G.v_in);

% Old doc from Vassilis

% [v_i, v_j, weights, ind_e] = GRAPH_ADJ2VEC(A):
%
% Convert adjacency matrix to edge vector form.
%
% For a graph with |V| nodes and |E| edges we have:
%
%       A: adjacency matrix                 |V| x |V|
%       v_i, v_j: indices of the edges      |E| x 1 
%
% that is:
% there is an edge from vertex v_i(e) to node v_j(e) for all e \in E
%
%
%INPUT:
% A: adjacency matrix with weights
%
%OUTPUTS:
%
% weights(e) = A( v_i(e), v_j(e) )
%            = A( ind_e(e) )
%
% weights will contain all edges only once, i.e. will use only the triu
% part of the adjacency matrix. It also handles the loops (diagonal
% elements of A).
%
%
%see also: graph_grad, graph_div
%
% code author: Vassilis Kalofolias
% date: August 2013













