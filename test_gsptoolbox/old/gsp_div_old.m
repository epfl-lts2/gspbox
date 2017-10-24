% div = GRAPH_DIV(S, A, ki, kj)
%
%Divergence of signal residing on edges of a graph. The output resides on
%nodes. The result should be such that 
%   graph_div(graph_grad(s)) = L * s,
%where L is the graph Laplacian associated with the adjacency matrix A.
%
%OUTPUT:
%   div: divergence of columns of S defined on vertices
%
%INPUTS:
%   S: signal defined on edges e
%   A: adjacency matrix (weighted or not)
%   ki, kj: indices of edges by graph_adj2vec
%
%
% TODO: does not handle normalized Laplacian case!! 
%
%
%note that the norm of the operator is
%
% ||D||^2_2 = ||L||_2, where G is the divergence operator and L is the
% Laplacian used.
%
%
%see also: graph_grad, graph_adj2vec, sgwt_laplacian
%
%code author: Vassilis Kalofolias
%date: Aug 2013

function div = gsp_div_old(G,s)


if ~strcmp(G.lap_type,'combinatorial') 
    error('Not implemented yet. However ask Nathanael it is very easy');
end


if size(s,1) ~= G.Ne
    error('Signal size not equal to number of edges');
end

k = size(s,2);

% What happens with loops? (diagonal elements of G, A)
%
% I think both div and grad is zero on the diagonal elements. But then is
% it consistent with the Laplacian? (it has to be div(grad(s)) = 2*L*s) for
% a signal s defined on the vertices.

div = zeros(G.N, k);

for ii = 1 : k
    tmp = sparse(G.v_in, G.v_out, s(:,ii), G.N,G.N);
    tmp = tril(tmp,1) - tril(tmp,1)';

    % the divergence 
    div(:,ii) = sum(tmp.* sqrt(G.W), 2);
end

end
