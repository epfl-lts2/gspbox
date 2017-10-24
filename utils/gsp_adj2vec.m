function [G] = gsp_adj2vec(G)
%GSP_ADJ2VEC Prepare the graph for the gradient computation
%   Usage: [G] = gsp_adj2vec(G)
%
%   Input parameters:
%       G   : Graph structure
%   Output parameters:
%       G   : Graph structure
%
%   This function converts adjacency matrix to edge vector form. It also
%   add the field G.Diff that is the sparse gradient matrix
%
%   See also: gsp_grad gsp_div

% Author: Nathanael Perraudin, Vassilis Kalofolias, Francesco Grassi
% Date  : 5 August 2014
% Testing: test_operators

if isfield(G,'Gm');
    G.Gm = gsp_adj2vec(G.Gm);
end

if G.directed
    % Decompose adjacency matrix W = Wsym + Wasym
    [Wsym,Wasym]=gsp_decompose_asymmatrix(G.W);
    
    % Find sym and asym edges
    [v_sym_i, v_sym_j, weights_sym] = find(Wsym);
    [v_asym_i, v_asym_j, weights_asym] = find(Wasym);
    
    
    G.v_out = v_asym_i;
    G.v_in = v_asym_j;
    G.weights = weights_asym;
    G.Ne = length(G.v_in);
    Diff_dir = gsp_grad_mat( G );
    
    G.v_out = v_sym_i;
    G.v_in = v_sym_j;
    G.weights = weights_sym;
    G.Ne = length(G.v_in);
    Diff_und = gsp_grad_mat( G );
    
    
    G.v_out = [v_sym_i;v_asym_i];
    G.v_in = [v_sym_j;v_asym_j];
    G.weights = [weights_sym;weights_asym];
    G.Ne = length(G.v_in);
    
    % Compute gradient matrix
    G.Diff = [Diff_dir;Diff_und];

    % Compute flow matrix
    G.Fo = [ Diff_dir;Diff_und  ];
    G.Fo(G.Fo>0) = 0; 
    
    G.Fi = [ Diff_dir;Diff_und  ];
    G.Fi(G.Fi<0) = 0;
    
    G.Adv = G.Diff'*G.Fo;
    
else
    % Keep each edge only once (they are duplicated!). Keep also loops.
    [v_i, v_j, weights] = find(tril(G.W));
    G.v_out = v_j;
    G.v_in = v_i;

    % the indices of the edges in the Adgacency matrix A:
    %G.ind_edges = sub2ind(size(G.W), G.v_in, G.v_out);
    G.weights = weights;         % |E| x 1
    G.Ne = length(G.v_in);

    G.Diff = gsp_grad_mat( G );
end

end




% Old code from Vasilis
% [v_i, v_j] = find(G.W);
% 
% % Keep each edge only once (they are duplicated!). Keep also loops.
% ind_keep = (v_i >= v_j);
% G.v_in = v_i(ind_keep);
% G.v_out = v_j(ind_keep);
% 
% % the indices of the edges in the Adgacency matrix A:
% G.ind_edges = sub2ind(size(G.W), G.v_in, G.v_out);
% G.weights = full(G.W(G.ind_edges));         % |E| x 1
% G.Ne = size(G.v_in);

% Old doc from Vassilis

% [v_i, v_j, weights, ind_e] = GRAPH_ADJ2VEC(A):
%
% Convert adjacency matrix to edge vector form.
%
% For a graph with |V| nodes and |E| edges we have:
%
%       A: adjacency matrix                 |V| x |V|
%       v_i, v_j: indices of the edges      |E| x 1 
%
% that is:
% there is an edge from vertex v_i(e) to node v_j(e) for all e \in E
%
%
%INPUT:
% A: adjacency matrix with weights
%
%OUTPUTS:
%
% weights(e) = A( v_i(e), v_j(e) )
%            = A( ind_e(e) )
%
% weights will contain all edges only once, i.e. will use only the triu
% part of the adjacency matrix. It also handles the loops (diagonal
% elements of A).
%
%
%see also: graph_grad, graph_div
%
% code author: Vassilis Kalofolias
% date: August 2013












