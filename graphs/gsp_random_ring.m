function [G]=gsp_random_ring(N)
%GSP_RANDOM_RING  Initialize a random ring graph
%   Usage:  G = gsp_random_ring(N);
%           G = gsp_random_ring();
%
%   Input parameters:
%         N     : Number of vertices. (default 64)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_ring(N)' 
%   weights are all equal to 1.
%
%   Example:
%
%          G = gsp_random_ring(64);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_random_ring.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Author : Nathanael Perraudin

if nargin < 1
   N = 64; 
end

G.N=N;

position = rand(N,1);
position = sort(position,'ascend');

weight = N*diff(position);
weightend = N * (1 + position(1) - position(end));

inds_j = (2:N)';
inds_i = (1:(N-1))';

G.W = sparse(inds_i, inds_j, weight, N,N);
G.W(N,1) = weightend;

G.W = G.W +G.W';


% Create coordinates
G.coords=[(cos(position*2*pi)),(sin(position*2*pi))];
G.limits=[-1,1,-1,1];


G.type = 'random-ring';


G = gsp_graph_default_parameters(G);

end

