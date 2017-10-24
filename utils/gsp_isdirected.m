function bool = gsp_isdirected(G)
%GSP_ISDIRECTED Check is the graph is directed
%   Usage: bool = gsp_isdirected(G);
%          bool = gsp_isdirected(W);
%
%   Input parameters
%       G       : Graph structure or square matrix
%   Output parameters
%       bool    : Boolean
%
%   This function test if the graph is directed. Alternatively, you can
%   give a square matrix and it tests if it is symetric. The function
%   returns 0 if the matrix is symetric and 1 otherwise!
%

% Author: Nathanael Perraudin
% Date  : 13 August 2014

if isstruct(G)
    W = G.W;
else
    W = G;
end

bool = sum(sum( abs(W - transpose(W))> eps(10) ))>0;

end
