function [ D ] = gsp_grad_mat( G )
%GSP_GRAD_MAT Gradient sparse matrix of the graph G
%   Usage:  D = gsp_gradient_mat(G);
%
%   Input parameters:
%       G   : Graph structure
%
%   Output parameters:
%       D   : Gradient sparse matrix
%
%   This function return the gradient matrix. To be more effiecient, call
%   the function:: 
%
%           G = gsp_adj2vec(G)
%
%   before this function.
%
%   Example:::
%
%       N = 40;
%       G = gsp_sensor(N);
%       G = gsp_adj2vec(G);
%       D = gsp_grad_mat(G);
%

% Author: Nathanael Perraudin
% Date  : 14 Mai 2014
% Testing: test_operators


if ~isfield(G,'v_in')
    G = gsp_adj2vec(G);
    warning(['GSP_GRADIENT_MAT: To be more efficient you should run: ',...
        'G = gsp_adj2vec(G); before using this proximal operator.']);
end

% if isfield(G,'Diff');
%     D = G.Diff;
%     return;
% end

% In case the graph has logical in the weights matrix
G.weights = double(G.weights);

if strcmp(G.lap_type,'combinatorial')
    n = G.Ne;
    Dr = [1:n 1:n];
    Dc(1:n) = G.v_in;
    Dc(n+1:2*n) = G.v_out;
    Dv(1:n) = sqrt(G.weights);
    Dv(n+1:2*n) = -sqrt(G.weights);
elseif strcmp(G.lap_type,'normalized')
    n = G.Ne;
    Dr = [1:n 1:n];
    Dc(1:n) = G.v_in;
    Dc(n+1:2*n) = G.v_out;
    Dv(1:n) = sqrt(G.weights./G.d(G.v_in));
    Dv(n+1:2*n) = -sqrt(G.weights./G.d(G.v_out));
%     [v_i, v_j, weights] = find(G.W);
%     n = length(v_i);
%     Dr = [1:n 1:n];
%     Dc(1:n) = v_i;
%     Dc(n+1:2*n) = v_j;
%     Dv(1:n) = sqrt(weights./G.d(v_i));
%     Dv(n+1:2*n) = -sqrt(weights./G.d(G.v_j));   
else
    error('Not implemented yet!')
end

    
D = sparse(Dr,Dc,Dv,n,G.N);

end

