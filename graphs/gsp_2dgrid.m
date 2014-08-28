function [G]=gsp_2dgrid(N,M)
%GSP_2dgrid  Initialize a 2 dimentional grid graph
%   Usage:  G=gsp_path(N);
%
%   Input parameters:
%         N     : Number of vertices along the first dimention (default 16)
%         M     : Number of vertices along the second dimention (default N)
%   Output parameters:
%         G     : Graph structure.
%
%   The 2d grid graph correspond the  graph used for the DCT2????
%
%   Example:
%
%          G = gsp_2dgrid(16);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_path
%
%   References:
%     G. Strang. The discrete cosine transform. SIAM review, 41(1):135-147,
%     1999.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_2dgrid.php

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

% Author: Nathanael Perraudin
% Date: 15 March 2014
% Testing: test_graph

if nargin <1
   N = 16; 
end

if nargin < 2
   M = N; 
end

G.N=N*M;

% Create weighted adjancency matrix
K = 2*(N-1);
J = 2*(M-1);
i_inds = zeros(K*M+J*N,1);
j_inds = zeros(K*M+J*N,1);
for ii = 1:M
    i_inds((ii-1)*K+(1:K)) = (ii-1)*N+[1:N-1,2:N]';
    j_inds((ii-1)*K+(1:K)) = (ii-1)*N+[2:N,1:N-1]';
end

for ii = 1:M-1
    i_inds(K*M+(ii-1)*2*N+(1:2*N)) = [((ii-1)*N+(1:N)),(ii*N+(1:N))]';
    j_inds(K*M+(ii-1)*2*N+(1:2*N)) = [(ii*N+(1:N)),((ii-1)*N+(1:N))]'; 
end


G.W = sparse(i_inds,j_inds,ones(K*M+J*N,1),N*M,N*M);

% Create coordinates
G.coords = [repmat((1:N)',M,1),...
        reshape(repmat((1:M),N,1),M*N,1) ];
G.plotting.limits = [0, N+1, 0, M+1];

G.type = '2d-grid';
G.plotting.vertex_size = 30;

G = gsp_graph_default_parameters(G);

end

