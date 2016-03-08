function G = gsp_stochastic_block_graph(N,k)
%GSP_STOCHASTIC_BLOCK_GRAPH  Create a stochastic block graph
%   Usage:  G = gsp_stochastic_block_graph( N );
%           G = gsp_stochastic_block_graph(N , k);
%
%   Input parameters:
%         N     : Number of nodes (default 5)
%         k     : Number of clusters (default 1024)
%   Output parameters:
%         G     : Graph structure.
%
%   Use the stochastic block model to create a graph.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_stochastic_block_graph.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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

% Author: Pierre Vandergheynst, Nathanael Perraudin
% Date  : 2 novemeber 2015



%% Stochastic Block Model generator
if nargin<1
    N = 1024; % number of nodes
end
if nargin<2
    k = 5; % number of clusters
end

%% partition with clusters of uniform random sizes
%  will generate clusters of similar sizes
z = randi(k,1,N); % partition vector
 
%% partition with clusters of homogenous sizes
%
% THIS PART IS NOT NECESSARY BECAUSE SLOW AND REPLACED
% BY THE SIMPLE "UNIFORM" CLUSTER ASSIGNMENT ABOVE
%
% L = ceil(N/k);
% for i=1:k-1
%     z((i-1)*L + 1:i*L) = i;
% end
% z((k-1)*L + 1:end) = k;
 
%% Link probability matrix
M = ones(k,k); % block probability matrix
p = 0.7; % proba link in cluster
q = 1-p; % proba link out of cluster
 
M = M * q/k;
M = M - diag(diag(M)) + diag(p*ones(1,k));
 
%% Generate adjacency
W = rand(N,N);
 
for i=1:N
    for j=1:N
        W(i,j) = ( W(i,j) <= M(z(i),z(j)) );
    end
end
 
 
 
[ G ] = gsp_graph_default_parameters( W );

end
