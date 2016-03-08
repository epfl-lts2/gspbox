function [ G ] = gsp_full_connected( N )
%GSP_FULL_CONNECTED  Create a fully connected graph
%   Usage:  G = gsp_full_connected(N);
%           G = gsp_full_connected();
%
%   Input parameters:
%         N     : Number of vertices (default 10)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_full_connected(N)' initializes a graph structure representing a
%   fully connected graph. All weight are set to 1. 
%
%   Example:
%
%          G = gsp_full_connected(5);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_full_connected.php

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



G.N = N;

W = ones(N);
W = W-diag(diag(W));

G.W = W;


% Create coordinates
G.coords = [(cos((0:N-1)*(2*pi)/N))',(sin((0:N-1)*(2*pi)/N))'];
G.plotting.limits = [-1,1,-1,1];

G.type = 'full';

G = gsp_graph_default_parameters(G);


end


