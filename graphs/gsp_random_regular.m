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
%   Example:::
%
%          G = gsp_random_regular(100,3)
%

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



