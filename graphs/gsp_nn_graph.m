function [ G ] = gsp_nn_graph(Xin, param)
%GSP_NN_GRAPH Create a nearest neighbors graph from a point cloud
%   Usage :  G = gsp_nn_graph( Xin );
%            G = gsp_nn_graph( Xin, param );
%
%   Input parameters:
%       Xin         : Input points
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_nn_graph( Xin, param )' creates a graph from positional data. The points are 
%   connected to their neighbors (either belonging to the k nearest 
%   neighbors or to the epsilon-closest neighbors. 
%
%   Example:::
%
%           P = gsp_pointcloud('bunny');
%           param.type = 'knn';
%           G = gsp_nn_graph(double(P), param);
%           gsp_plot_graph(G);
%
%   Additional parameters
%   ---------------------
%
%   * *param.type*      : ['knn', 'radius']   the type of graph (default 'knn')
%   * *param.use_flann* : [0, 1]              use the FLANN library (default 0)
%   * *param.use_full*  : [0, 1] - Compute the full distance matrix and then
%     sparsify it (default 0) 
%   * *param.center*    : [0, 1]              center the data (default 0)
%   * *param.rescale*   : [0, 1]              rescale the data on a 1-ball (def 0)
%   * *param.sigma*     : float               the variance of the distance kernel
%   * *param.k*         : int                 number of neighbors for knn (def 10)
%   * *param.epsilon*   : float               the radius for the range search
%   * *param.use_l1*    : [0, 1]              use the l1 distance (def 1)
%   * *param.symmetrize_type*: ['average','full'] symmetrization type (default 'full')
%   * *param.min_weight*: float               constant additive weight for each edge (default 0) 
%   * *param.zero_diagonal*: [0, 1]           zero out the diagonal (default 1)
%   * *param.weight_kernel*: function         edge-weighting kernel (default gaussian)
%
%
%   See also: gsp_nn_distanz gsp_pointcloud
%

% Author: Johan Paratte, Nathanael Perraudin
% Date: 16 June 2014
% Testing: test_rmse

    if nargin < 2
    % Define parameters
        param = {};
    end
    
    %Parameters
    if ~isfield(param, 'type'), param.type = 'knn'; end
    if ~isfield(param, 'use_flann'), param.use_flann = 0; end
    if ~isfield(param, 'center'), param.center = 0; end
    if ~isfield(param, 'zero_diagonal'), param.zero_diagonal = 1; end
    if ~isfield(param, 'rescale'), param.rescale = 0; end
    if ~isfield(param, 'k'), param.k = 10; end
    if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
    if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
    if ~isfield(param, 'target_degree'), param.target_degree = 0; end;
    if ~isfield(param, 'symmetrize_type'), param.symmetrize_type = 'average'; end
    if ~isfield(param, 'min_weight'), param.min_weight = 0; end
    if ~isfield(param, 'light'); param.light = 0; end
    if ~isfield(param, 'weight_kernel'); param.weight_kernel = @(x, sigma) exp(-x/sigma); end
    
    paramnn = param;
    paramnn.k = param.k +1;
    kdist = @(x1,x2) gsp_nn_distanz(x1',x2',paramnn);
    [indx, indy, dist, Xout, ~, epsilon] = kdist(Xin,Xin);
    Xout = transpose(Xout);
    switch param.type
        case 'knn'
            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = mean(dist); end
            else
                if ~isfield(param, 'sigma'), param.sigma = mean(dist)^2; end
            end
        case 'radius'
            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = epsilon/2; end
            else
                if ~isfield(param, 'sigma'), param.sigma = epsilon.^2/2; end
            end
        otherwise
            error('Unknown graph type')
    end
    
    
    if param.use_l1
        Wmat = @(indx,indy,dist,n,m) sparse(indx, indy, param.min_weight+double(param.weight_kernel(dist, param.sigma)), n, m);
    else
        Wmat = @(indx,indy,dist,n,m) sparse(indx, indy, param.min_weight+double(param.weight_kernel(dist.^2, param.sigma)), n, m);
    end
    
    n = size(Xin,1);
    W = Wmat(indx,indy,dist,n,n);
    
    Wkernel = @(x) create_kernel(Xin,x,kdist,Wmat);
    
    if param.zero_diagonal
        % We need zero diagonal
        W(1:(n+1):end) = 0;     % W = W-diag(diag(W));
    end
    
    % Computes the average degree when using the epsilon-based neighborhood
    if (strcmp(param.type,'radius'))
        text = sprintf('Average number of connection = %d', nnz(W)/size(W, 1));
        disp(text);
    end

    % Sanity check
    if size(W,1) ~= size(W,2), error('Weight matrix W is not square'); end
    
    % Symmetry checks
    %if issymmetric(W)
    if (norm(W - W', 'fro') == 0)
        disp('The matrix W is symmetric');
    else
         W = gsp_symmetrize(W,param.symmetrize_type);
    end
    
    %Fill in the graph structure
    G.N = n;
    G.W = W;
    G.coords = Xout;
    G.Wkernel = Wkernel;
    %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
    if param.use_l1
        G.type = 'nearest neighbors l1';
    else
        G.type = 'nearest neighbors';
    end
    %G.vertex_size=30;
    G.sigma = param.sigma;
    if param.light
        G = gsp_graph_lightweight_parameters(G);
    else
        G = gsp_graph_default_parameters(G);
    end
end

function W = create_kernel(Xin,x,kdist,Wmat)
[indx, indy, dist] = kdist(Xin,x);
n = size(Xin,1);
m = size(x,1);
W = Wmat(indx,indy,dist,n,m)';

end
