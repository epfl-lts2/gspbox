function [G]=gsp_logo()
%GSP_LOGO Initialize a graph with the GSP logo
%   Usage:  G = gsp_logo();
%
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_logo()' initializes a graph structure containing the GSP logo
%
%   Example:::
%
%          G = gsp_logo();
%          gsp_plot_graph(G);
%

% Author : Johan Paratte
% Test: 

load logogsp

G.W = W;
G.coords = coords;
G.info.idx_g = idx_g;
G.info.idx_s = idx_s;
G.info.idx_p = idx_p;
%G.plotting.vertex_color = [200 136/255.0 204/255.0];
%G.plotting.edge_color = [0 136/255.0 204/255.0];
G.plotting.vertex_size = 20;
G.limits = [0 640 -400 0];
G = gsp_graph_default_parameters(G);

end


