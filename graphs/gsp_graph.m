function [G] = gsp_graph(W,coords,limits)
%GSP_GRAPH  Initialize a graph from a weight matrix
%   Usage:  G = gsp_graph(N);
%
%   Input parameters:
%         W     : Weight matrix
%         coords: Coordonates of the points (optional)
%         limits: limits for the coordonates (optional)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_graph(W,coords,limits)' initializes a graph structure with W as
%   weight matrix.
%
%   Example:
%
%          W = rand(10);
%          W = W - diag(diag(W));
%          W = (W + W')/2;
%          G = gsp_graph(W);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_graph.php

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
% Date: 16 March 2014

gsp_check_weights(W);

G.W = W;

% Create coordinates
if nargin > 1
    G.coords = coords;
end
if nargin > 2
    G.plotting.limits = limits;
end

G.type = 'from weight';

G = gsp_graph_default_parameters(G);

end

