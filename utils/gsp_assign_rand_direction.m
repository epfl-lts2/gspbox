function G = gsp_assign_rand_direction(G,p)
% GSP_ASSIGN_RAND_DIRECTION Assign random direction to p percent edges of the undirected graph G
%
%   Usage: G = gsp_assign_rand_direction(G,p)
%
%   Input parameters:
%       G   : Undirected Graph structure
%       p   : Percentage of edges
%   Output parameters:
%       G   : Directed Graph structure
%
%   Assign random direction to p percent edges of the undirected graph G
%

% Author: Francesco Grassi
% Date  : 4 July 2014
% Testing: test_operators

if nargin<2
    p = 1;
end
W = G.W;

[i]=find(triu(W));

k = randperm(length(i),round(p*length(i)));

W(i(k))=0;

G = gsp_update_weights(G,W);

end