function coords = gsp_laplacian_eigenmaps(G, dim, param)
%GSP_LAPLACIAN_EIGENMAPS Laplacian eigenmaps embedding 
%   Usage: coords = gsp_laplacian_eigenmaps(G, dim);
%          coords = gsp_laplacian_eigenmaps(G, dim, param);
%
%   Input parameters
%         G      : Graph
%         dim    : Dimension of the embedding
%         param  : Structure of optional paramters
%   Output parameters
%         coords : Coordinates of the embedding
% 
%   This function uses the weight matrix of a graph G, in order to compute
%   a *dim* -dimensional embedding (output coordinates). The algorithm used
%   is Laplacian eigenmaps. Warning, this function might not work if the
%   graph is not connected.
%
%   *param* is a structure optional parameters:
%     
%   * *param.tol* : Tolerance for the spectral gap (default 1e-6).   
%
%
%   See also: gsp_update_coordinates gsp_isomap gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%
%   References: belkin2001laplacian

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015


if nargin<3
    param = struct;
end

if ~isfield(param,'tol'), param.tol = 1e-6; end


% Laplacian Eigenmaps
options.disp = 0;
options.isreal = 1;
options.issym = 1;
% only need bottom (no_dims + 1) eigenvectors
if strcmp(G.type,'normalized') % Bug G.type
    [coords, e] = eigs(G.L, dim + 1, param.tol, options);
else
    [coords, e] = eigs(G.L, sparse(diag(G.d)), dim + 1, param.tol, options);
end
e = diag(e);

if nnz(e < param.tol) > 1
    warning('Multiple zero eigenvalues')
end

[~, ind] = sort(e, 'ascend');
coords = coords(:, ind);
coords = coords(:, 2:end);
end

