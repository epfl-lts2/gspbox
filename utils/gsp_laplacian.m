function L = gsp_laplacian(W, lap_type)
%GSP_LAPLACIAN compute the graph Laplacian from undirected adjacency matrix
%   Usage: L = gsp_laplacian(W, laplacian_type);
%          G = gsp_laplacian(W);
%
%   Input parameters:
%       W   : Weighted adjacency matrix
%       laplacian_type: Type of laplacian: 'combinatorial' or 'normalized'
%   Output parameters:
%       L   : Graph Laplacian
%
%   This function creates the graph laplacian given a weighted adjacency
%   matrix.
%
%   The variable *laplacian_type* contains the different laplacian types.
%   Available laplacian types:
%
%   * *combinatorial*: Non normalized laplacian (default).
%
%     ..   L =  D  - W 
%
%     .. math:: L= D - W
%   * *normalized*: Normalized laplacian
%     .. L_n = I - D^-0.5 W D^-0.5 
%
%     .. math:: L_{n}=I - D^{-\frac{1}{2}}WD^{-\frac{1}{2}}
%   * *none*: No laplacian (used for compatibility reasons).
%
%   For directed graphs, see gsp_create_laplacian.
%
%   References: chung2005laplacians
%
% see also: gsp_create_laplacian

% Author: Vassilis Kalofolias
% Date  : June 2016


if nargin < 2
    lap_type = 'combinatorial';
end

n = size(W,1);

d = sum(W,2);
switch lap_type
    case 'combinatorial'
        L = diag(d)-W;
    case 'normalized'
        if issparse(W)
            Dn = diag(d.^(-0.5));
            % CAREFUL: a zero row/column in W is followed by an inf in d,
            % but their multiplication gives zero again (because of using
            % sparse).
            L = speye(n) - Dn * W * Dn;
        else
            ind = d>0;
            L = - W;
            % Ln = D^(-1/2) L D^(-1/2)
            L(ind, :) = bsxfun(@times, L(ind, :), 1./sqrt(d(ind)));
            L(:, ind) = bsxfun(@times, L(:, ind), 1./sqrt(d(ind))');
            % put back diagonal to identity
            % Note: for disconnected nodes we should still have 1 on diagonal
            % (limit of L for W -> 0)
            L(1:n+1:end) = 1;
    %         L(ind, ind) = speye(nnz(ind))-Dn*W(ind, ind)*Dn;
        end
    case 'randomwalk'
        L = speye(n) - diag(d.^(-1))*W;

    case 'none'
        L = sparse(0);
    otherwise
        error(' Unknown laplacian type')
end
