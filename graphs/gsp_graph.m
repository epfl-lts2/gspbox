function G = gsp_graph(W, coords, limits)
%GSP_GRAPH  Create a graph given weighted adjacency matrix
%   Usage:  G = gsp_graph(W);
%           G = gsp_graph(W, coords);
%           G = gsp_graph(W, coords, limits);  
%
%   Input parameters:
%         W     : (n by n) Weighted adjacency matrix
%         coords: (n by 2) or (n by 3) Coordinates of the points (optional)
%         limits: limits for the coordinates (optional)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_graph(W, coords, limits)' initializes a graph structure with W as
%   weight matrix.
%
%   Example:::
%
%          W = rand(10);
%          W = W - diag(diag(W));
%          W = (W + W')/2;
%          G = gsp_graph(W);
%



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
