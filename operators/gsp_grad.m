function gr = gsp_grad(G,s)
%GSP_GRAD Graph gradient
%   Usage: gr = gsp_grad(G,s)
%   
%   Input parameters:
%       G   : Graph structure
%       s   : Signal living on the nodes
%   Output parameters:
%       gr  : Gradient living on the edges
%
%   For the non normalized Laplacian, the gradient of the node signal f
%   evaluated at the edge linking x and y is given by:
%
%           grad f (x,y) = sqrt( w(x,y) ) ( f(x) - f(y) )
%
%   Before using this function, you need to call the function: 
%
%           G = gsp_adj2vec(G)
%
%   See also: gsp_div gsp_adj2vec
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_grad.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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


if ~strcmp(G.lap_type,'combinatorial') 
    error('Not implemented yet. However ask Nathanael it is very easy');
end



D = gsp_grad_mat(G);

gr = D*double(s);

if isa(s,'single')
   gr = single(di); 
end



end



%% Old vassilis code and comment
% % the gradient is defined from the nodes on the edges 
% %grad = (s(ki) - s(kj)) .* (sqrt(w_e));
% G = bsxfun(@times, (s(G.v_in, :) - s(G.v_out, :)), sqrt(G.weights));      % |E| x k




% G = GRAPH_GRAD(S, ki, kj, w_e)
%
%Gradient of signal residing on nodes of a graph. The output resides on
%edges. The result should be such that
%   graph_div(graph_grad(s)) = L * s,
%where L is the graph Laplacian associated with the adjacency matrix A.
%
%INPUTS:
%   S:  signal defined on nodes v
%   ki, kj, w_e: indices and weights of edges (computed by graph_adj2vec)
%
%
%OUTPUT:
%   G:  gradient of columns of S (g defined on edges, has the same size as
%       inputs ki, kj and w_e) 
%
%
% TODO: does not handle normalized Laplacian case!! 
%
%
%
% note that the norm of the operator is
%
% ||G||^2_2 = ||L||_2, where G is the gradient operator and L is the
% Laplacian used.
%
%
%see also: graph_div, graph_adj2vec, sgwt_laplacian
%
%code author: Vassilis Kalofolias
%date: Aug 2013

% Note that the sign of the gradient above depends on the "direction" of
% the edge that is arbitrary here. We have computed the gradient only on
% the edges counted once, while in the adjacency matrix we had them
% duplicated (undirected graph). As a final output we could have a matrix
% in the size of A that has a zero diagonal and is antisymmetric (!).

% the form above is however more compact and we will keep it.





























