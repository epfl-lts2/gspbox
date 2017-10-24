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
%   Example:::
%
%          G = gsp_path(16);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_graph
%
%   References: strang1999discrete

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
