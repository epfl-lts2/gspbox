function G = gsp_path(N)
%GSP_PATH  Initialize a path graph
%   Usage:  G = gsp_path(N);
%           G = gsp_path();
%
%   Input parameters:
%         N     : Number of vertices (default 32).
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_path(N)' initializes a graph structure of a path graph.
%
%   The path graph correspond the  graph used for the DCT. See references
%   for more informations.
%
%   Example:
%
%          G = gsp_path(16);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_graph
%
%   References:
%     G. Strang. The discrete cosine transform. SIAM review, 41(1):135--147,
%     1999.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_path.php

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

% Author David I Shuman, Nathanael Perraudin
% Date: 15 March 2014

if nargin < 1
    N = 16;
end


% Create weighted adjancency matrix
i_inds = [1:N-1, 2:N];
j_inds = [2:N, 1:N-1];
G.W = sparse(i_inds, j_inds, ones(1,2*(N-1)), N, N);

% Create coordinates
G.coords = [(1:N)', zeros(N,1)];
G.plotting.limits = [0, N+1, -1, 1];

G.type = 'path';

G = gsp_graph_default_parameters(G);

end

