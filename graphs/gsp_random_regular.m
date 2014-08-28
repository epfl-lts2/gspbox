function [ G ] = gsp_random_regular( N,k )
%GSP_RANDOM_REGULAR Create a random regualar graph
%   Usage:  G = gsp_random_regular( N,k )
%           G = gsp_random_regular( N )
%           G = gsp_random_regular();
%
%   Input parameters:
%         N     : Number of nodes (default 64)
%         k     : Number of connection of each nodes (default 6)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_random_regular( N,k )' initializes create a graph structure
%   containing the weighted adjacency matrix (G.W), the number of vertices
%   (G.N) for a random regular graph. All edge weights are equal to 1. 
% 
%   The random regular graph has the property that every nodes is connected
%   to 'k' other nodes.
%
%   Example:
%
%          G = gsp_random_regular(100,3)
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_random_regular.php

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

% Author : Nathanael Perraudin
 

% Optional input arguments
if nargin < 2 
    k = 6;
end

if nargin < 1
   N = 64; 
end

G.type = 'random_regular';

G.W = createRandRegGraph(N,k); % Use downloaded code 'createRandRegGraph' 

G = gsp_graph_default_parameters(G);

end




