function coords = gsp_isomap(G, dim, param)
%GSP_ISOMAP isomap
%   Usage: coords = gsp_isomap(G, dim, param);
%          coords = gsp_isomap(G, dim);
%
%   Input parameters
%         G         : Graph
%         dim       : Dimensionality of the embedding
%         param     : Structure of optional parameters
%
%   Output parameters
%         coords    : Coordinates of the embedding
%
%   This function uses the weight matrix of a graph G, in order to compute
%   a *dim* -dimensional embedding (output coordinates). The algorithm used
%   is Isomap. Warning, this function might not work if the
%   graph is not connected.
%
%   *param* is a structure with optional parameters:
%
%   * *param.tol*    : Tolerance for the spectral gap (default 1e-6).
%   * *param.kernel* : The kernel used to create the graph weight matrix:
%     + 'exp'        : exponential kernel ($e^(frac{-d_{ij}}{sigma^2})$)
%     + '1/x'        : inverse of x kernel ($frac{1}{sigma+d_{ij}}$)
%     + '1/x^2'      : inverse of x^2 kernel ($frac{1}{(sigma+d_{ij})^2}$)
%     + 'resistance' : Resistance distance.
%   * *param.k*      : Max number of nearest neighbors. If not defined, the
%
%   References: tenenbaum2000global
%
%   See also: gsp_update_coordinates gsp_laplacian_eigenmaps gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%
%   References: tenenbaum2000global
%

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

if nargin<3
    param = struct;
end


if ~isfield(param,'kernel'), param.kernel = 'exp'; end

D = gsp_weight2distance(G,param.kernel);

% Shortest paths
d = zeros(G.N);
for ii = 1:G.N
    d(ii, :) = dijkstra(D, ii); 
end
D = (d.^2);



if any(isinf(D))
    warning('Disconnected graph!');
    D(isinf(D)) = max(D(not(isinf(D)))) * 10;
end

% Perform MDS on Dijkstra Distance matrix D
B = 0.5*(repmat(sum(D,2),1,length(D))/length(D)...
    + repmat(sum(D,2).',length(D),1)/length(D) ...
    - D - repmat(sum(sum(D,1)),...
    length(D),length(D))/length(D).^2);


[coords, e] = eigs(B);
e = diag(e);
[e , ind] = sort(e, 'descend');
coords = coords(:,ind(1:dim));
coords = coords * diag(sqrt(e(1:dim)));

end