function [ G ] = gsp_update_weights( G, W )
%GSP_UPDATE_WEIGHTS update weight matrix of the graph G
%   Usage: G = gsp_update_weight( G, W );
%
%   Input parameters
%       G   : Graph
%       W   : W new weight matrix
%   Output parameters
%       G   : Graph
% 
%   This function will update the weight of graph and recompute the
%   Laplacian the new Laplacian. It will also recompute the value G.lmax*
%   if it is present in the graph.
%
%   For now this function does not recompute the Fourier basis!
%
%   Example:
%          
%          G.W = Wnew;
%          G = gsp_graph_default_parameters( G );
%
%          
%   See also: gsp_graph_default_parameters gsp_graph
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_update_weights.php

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

% Author: Nathanael Perraudin
% Date  : 19.11.2014




G.W = W;

if isfield(G,'U')
    G = rmfield(G,'U');
end

if isfield(G,'e')
    G =rmfield(G,'e');
end

if isfield(G,'mu')
    G =rmfield(G,'mu');
end

G = gsp_graph_default_parameters(G);

if isfield(G,'lmax')
    G = rmfield(G,'lmax');
    G = gsp_estimate_lmax(G);
end


end


