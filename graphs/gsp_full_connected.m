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
%   Example:::
%
%          G = gsp_full_connected(5);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring



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

