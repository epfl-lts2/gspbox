function G = gsp_rmse_mv_graph(X,param)
%GSP_RMSE_MV_GRAPH Root meean square error missing value graph
%   Usage :  G = gsp_rmse_mv_graph( Xin );
%            G = gsp_rmse_mv_graph( Xin, param );
%
%   Input parameters:
%       Xin         : Input points (with missing value (nan))
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_rmse_mv_graph(X,param)' creates a graph from positional data. The points are 
%   connected to their neighbors (either belonging to the k nearest 
%   neighbors or to the epsilon-closest neighbors. This function ignore all
%   nan value. But it is much slower than |gsp_nn_graph|.
%
%   Additional parameters
%   ---------------------
%
%   * *param.type*      : ['knn', 'radius']   the type of graph (default 'knn')
%   * *param.sigma*     : float               the variance of the distance kernel
%   * *param.k*         : int                 number of neighbors for knn
%   * *param.epsilon*   : float               the radius for the range search
%   * *param.symmetrize_type*: ['average','full'] symmetrization type (default 'full')
%   * *param.center*    : [0, 1]              center the data
%   * *param.rescale*   : [0, 1]              rescale the data (in a 1-ball)

% Author: Nathanael Perraudin
% Date : 12 March 2015
% Testing: test_rmse

if nargin<2
    param = struct;
end

if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'k'), param.k = 10; end
if ~isfield(param, 'symmetrize_type'), param.symmetrize_type = 'average'; end
if ~isfield(param, 'type'), param.type = 'knn'; end
if ~isfield(param, 'center'), param.center = 1; end
if ~isfield(param, 'rescale'), param.rescale = 1; end


[N, d] = size(X);


if N> 5000
    error('This code is not optimized and cannot handle such a big graphs!');
end

Xout = X;
%Center the point cloud
if param.center
    for ii = 1:d
        mask = 1-isnan(X(:,ii));
        mask = mask .* (1:N)';
        mask = find(mask);
        Xout(mask,ii) = X(mask,ii) - mean(X(mask,ii));
    end
end

%Rescale the point cloud
if param.rescale
    bounding_radius = 0;
    for ii = 1:d
        mask = 1-isnan(X(:,ii));
        mask = mask .* (1:N)';
        mask = find(mask);
        bounding_radius = bounding_radius+ abs(max(Xout(mask,ii)) - min(Xout(mask,ii)))^2;
    end
    bounding_radius = 0.5 * sqrt(bounding_radius);
    scale = nthroot(N, min(d, 3)) / 10;
    Xout = Xout .* (scale / bounding_radius);
end



p.verbose = param.verbose;
C = gsp_rmse_mv(transpose(Xout),p);
C = C*sqrt(d);



% sparsification
switch param.type
    case 'knn'
        for ii=1:size(C,1)
            [~,sortIndex] = sort(C(ii,:),'ascend');  %# Sort the values in
                                                     %#   descending order
            C(ii,sortIndex(param.k+2:end)) = 0;
        end
        if ~isfield(param, 'sigma'), param.sigma = mean(C(C>0))^2; end
        
        W = spfun(@(x) exp(-x.^2/param.sigma),sparse(C));
        W = W-diag(diag(W));

        W = gsp_symmetrize(W,param.symmetrize_type);
    case 'radius'
        if ~isfield(param, 'sigma'), param.sigma = param.epsilon.^2/2; end

        C(C>param.epsilon) = 0;
        W = spfun(@(x) exp(-x.^2/param.sigma),sparse(C));
        W = W-diag(diag(W));
    otherwise
       error('Unknown type : allowed values are knn, radius');
end


% Create the graph
G = gsp_graph(W,Xout);
G.type = 'RMSE MV';
G.sigma = param.sigma;


        
end
