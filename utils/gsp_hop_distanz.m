function d = gsp_hop_distanz(G,i,j)
%GSP_HOP_DISTANZ Compute the hop distance between two node
%   Usage:  d = gsp_hop_distanz(G,i,j);
%
%   Input parameters:
%       G   : Graph
%       i   : node
%       j   : node
%   Output parameters:
%       d   : hop distanz
%
%   This code computes the hop distance between node i and node j. It uses
%   a naive greedy algorithm and has to be improved.
%

% Author: Nathanael Perraudin
% Date  : 15 septembre 2015
% Testing: test_gsp_hope_distanz

M = double(logical(G.W));
s = zeros(G.N,1);
s(i) = 1;
s = logical(s);
d = 0;
while s(j)==0 && d<=G.N+1
    d = d+1;
    s = logical(M*double(s)) +s;
end

if d == G.N+1
    d = inf;
end


end