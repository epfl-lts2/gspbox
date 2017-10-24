function [G]=gsp_random_ring(N)
%GSP_RANDOM_RING  Initialize a random ring graph
%   Usage:  G = gsp_random_ring(N);
%           G = gsp_random_ring();
%
%   Input parameters:
%         N     : Number of vertices. (default 64)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_ring(N)' 
%   weights are all equal to 1.
%
%   Example:::
%
%          G = gsp_random_ring(64);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%

% Author : Nathanael Perraudin

if nargin < 1
   N = 64; 
end

G.N=N;

position = rand(N,1);
position = sort(position,'ascend');

weight = N*diff(position);
weightend = N * (1 + position(1) - position(end));

inds_j = (2:N)';
inds_i = (1:(N-1))';

G.W = sparse(inds_i, inds_j, weight, N,N);
G.W(N,1) = weightend;

G.W = G.W +G.W';


% Create coordinates
G.coords=[(cos(position*2*pi)),(sin(position*2*pi))];
G.limits=[-1,1,-1,1];


G.type = 'random-ring';


G = gsp_graph_default_parameters(G);

end
