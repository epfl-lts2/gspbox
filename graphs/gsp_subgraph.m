function [ subG ] = gsp_subgraph( G,c )
%GSP_SUBGRAPH Create a subgraph from G
%   Usage: [ subG ] = gsp_subgraph( G,c )
%
%   Input parameters:
%       G   : Original graph
%       c   : Node to keep
%
%   Output parameters:
%       subG: Subgraph
%
%   This function create a subgraph from G taking only the node in c.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_subgraph.php

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
% Date : Mai 2015


subG.W = G.W(c,c);
subG.N = length(c);
if isfield(G,'type')
    subG.type = ['sub-',G.type];
end

if isfield(G,'directed')
    subG.directed = G.directed;
end

if isfield(G,'lap_type')
    subG.lap_type = G.lap_type;
end

if isfield(G,'coords')
    subG.coords = G.coords(c,:);
end

if isfield(G,'limits')
    subG.limits = G.limits;
end


if isfield(G,'edge_width')
    subG.edge_width = G.edge_width;
end


if isfield(G,'edge_color')
    subG.edge_color = G.edge_color;
end


if isfield(G,'edge_style')
    subG.edge_style = G.edge_style;
end


if isfield(G,'vertex_size')
    subG.vertex_size = G.vertex_size;
end


if isfield(G,'vertex_color')
    subG.vertex_color = G.vertex_color;
end


if isfield(G,'vertex_edge_color')
    subG.vertex_edge_color = G.vertex_edge_color;
end

subG = gsp_graph_default_parameters(subG);

end


