function [ G ] = gsp_nn_hypergraph( Xin, param )
%GSP_NN_HYPERGRAPH Create a nearest neighbors hypergraph from a point cloud
%   Usage :  G = gsp_nn_hypergraph( Xin );
%            G = gsp_nn_hypergraph( Xin, param );
%
%   Input parameters:
%       Xin         : Input points
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   Example:::
%
%           P = rand(100,2);
%           G = gsp_nn_hypergraph(P)
%
%   Additional parameters
%   ---------------------
%
%   * *param.use_flann* : [0, 1]              use the FLANN library
%   * *param.center*    : [0, 1]              center the data
%   * *param.rescale*   : [0, 1]              rescale the data (in a 1-ball)
%   * *param.sigma*     : float               the variance of the distance kernel
%   * *param.k*         : int                 number of neighbors for knn
%
%   See also: gsp_nn_graph
%

% Author: Nathanael Perraudin
% Date: 21 October 2015
% Testing: test_rmse

%   * *param.type*      : ['knn', 'radius']   the type of graph (default 'knn')
%   * *param.epsilon*   : float               the radius for the range search
%   * *param.use_l1*    : [0, 1]              use the l1 distance

    if nargin < 2
    % Define parameters
        param = {};
    end
    
    %Parameters
%     if ~isfield(param, 'type'), param.type = 'knn'; end
    if ~isfield(param, 'use_flann'), param.use_flann = 0; end
    if ~isfield(param, 'center'), param.center = 0; end
    if ~isfield(param, 'rescale'), param.rescale = 1; end
    if ~isfield(param, 'k'), param.k = 10; end

    param.type = 'knn';
%     if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
%     if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
%     if ~isfield(param, 'target_degree'), param.target_degree = 0; end;
    paramnn = param;
%     paramnn.k = param.k +1;
    [indx, ~, dist] = gsp_nn_distanz(Xin',Xin',paramnn);
    
%     switch param.type
%         case 'knn'
%             if param.use_l1
%                 if ~isfield(param, 'sigma'), param.sigma = mean(dist); end
%             else
                if ~isfield(param, 'sigma'), param.sigma = mean(dist)^2; end
%             end
%         case 'radius'
%             if param.use_l1
%                 if ~isfield(param, 'sigma'), param.sigma = epsilon/2; end
%             else
%                 if ~isfield(param, 'sigma'), param.sigma = epsilon.^2/2; end
%             end
%         otherwise
%             error('Unknown graph type')
%     end
    

    w = exp(-dist.^2/param.sigma);
    G.N = size(Xin,1);
    G.Ne = G.N;
    G.W = sparse(G.N,G.Ne);
    G.E = cell(G.Ne,1);
    k = param.k;
    for ii = 1:G.Ne
        edge = indx((1:k)+(ii-1)*k);
        G.E{ii} = edge;
        % Here we use H for HW...
        G.W(edge,ii) = sqrt(sum(w(edge)));
    end
    G.hypergraph = 1;
    G.directed = 0;
    
    %Fill in the graph structure
    G.coords = Xin;

    G.type = 'Nearest neighboors hypergraph';
    G.lap_type = 'normalized';
    G.sigma = param.sigma;
    G = gsp_graph_default_parameters(G);
end

