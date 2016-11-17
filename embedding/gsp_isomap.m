function coords = gsp_isomap(G, dim, param)
%GSP_ISOMAP isomap
%   Usage: coords = gsp_isomap(G, dim, param);
%          coords = gsp_isomap(G, dim);
%
%   Input parameters
%         G         : Graph
%         dim       : Dimensionality of the embedding
%         param     : Structure of optional parameters
%
%   Output parameters
%         coords    : Coordinates of the embedding
%
%   This function uses the weight matrix of a graph G, in order to compute
%   a dim -dimensional embedding (output coordinates). The algorithm used
%   is Isomap. Warning, this function might not work if the
%   graph is not connected.
%
%   param is a structure with optional parameters:
%
%    param.tol    : Tolerance for the spectral gap (default 1e-6).
%    param.kernel : The kernel used to create the graph weight matrix:
%     + 'exp'        : exponential kernel (e^(frac{-d_{ij}}{sigma^2}))
%     + '1/x'        : inverse of x kernel (frac{1}{sigma+d_{ij}})
%     + '1/x^2'      : inverse of x^2 kernel (frac{1}{(sigma+d_{ij})^2})
%     + 'resistance' : Resistance distance.
%    param.k      : Max number of nearest neighbors. If not defined, the
%
%   References:
%     J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric
%     framework for nonlinear dimensionality reduction. Science,
%     290(5500):2319--2323, 2000.
%     
%
%   See also: gsp_update_coordinates gsp_laplacian_eigenmaps gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%
%   References:
%     J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric
%     framework for nonlinear dimensionality reduction. Science,
%     290(5500):2319--2323, 2000.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/embedding/gsp_isomap.php

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

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

if nargin<3
    param = struct;
end


if ~isfield(param,'kernel'), param.kernel = 'exp'; end

D = gsp_weight2distance(G,param.kernel);

% Shortest paths
d = zeros(G.N);
for ii = 1:G.N
    d(ii, :) = dijkstra(D, ii); 
end
D = (d.^2);



if any(isinf(D))
    warning('Disconnected graph!');
    D(isinf(D)) = max(D(not(isinf(D)))) * 10;
end

% Perform MDS on Dijkstra Distance matrix D
B = 0.5*(repmat(sum(D,2),1,length(D))/length(D)...
    + repmat(sum(D,2).',length(D),1)/length(D) ...
    - D - repmat(sum(sum(D,1)),...
    length(D),length(D))/length(D).^2);


[coords, e] = eigs(B);
e = diag(e);
[e , ind] = sort(e, 'descend');
coords = coords(:,ind(1:dim));
coords = coords * diag(sqrt(e(1:dim)));

end
