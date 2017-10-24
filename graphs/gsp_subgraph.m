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

