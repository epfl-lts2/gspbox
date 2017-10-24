function [ G ] = gsp_knn_classify_graph( Xl, Xu, param )
%GSP_KNN_CLASSIFY_GRAPH Create a weighted knn classifier graph
%   Usage :  G = gsp_knn_classify_graph( Xl, Xu );
%            G = gsp_knn_classify_graph( Xl, Xu, param );
%
%   Input parameters:
%       Xl         : Labeled points
%       Xu         : Unlabeled points
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   This function create a special graph from labeled and unlabeled data.
%   Solving a tikonow regression problem on this graph is equivalent to
%   perform a knn regression algorithm. This function is done for
%   comparison between methods.
%
%   Additional parameters
%   ---------------------
%
%   * *param.use_flann* : [0, 1]    use the FLANN library
%   * *param.k*         : int       number of neighbors for knn
%   * *param.use_l1*    : [0, 1]    use the l1 distance
%   * *param.weighted*  : [0, 1]    use a weighted graph
%   * *param.sigma*     : double    
%
%   See also: gsp_nn_graph gsp_knn_classify_weighted_graph
%

% Author: Nathanael Perraudin
% Date: 24 June 2015
% Testing: test_graph_ml

    if nargin < 3
    % Define parameters
        param = {};
    end
    
    %Parameters
    if ~isfield(param, 'use_flann'), param.use_flann = 0; end
    if ~isfield(param, 'k'), param.k = 10; end
    if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
    if ~isfield(param, 'weighted'), param.weighted = 0; end

    % test if the binaries of flann are working
    if param.use_flann
       try
            paramsflann.algorithm = 'kdtree';
            paramsflann.checks = 32;
            paramsflann.trees = 1;
            tmp = rand(100,10);
            [NN] = flann_search(tmp', tmp', 3, paramsflann); %#ok<NASGU>
       catch
            warning('Flann not compiled, going for the slow algorithm!')
            param.use_flann = 0;
       end
    end
    
    
    k = param.k;
    
    N = size(Xu,1);
    M = size(Xl,1);
    
    spi = zeros(N*k,1);
    spj = zeros(N*k,1);
    spv = zeros(N*k,1);

    %Find kNN for each point in X (Using a kdtree)
    if param.use_flann
        if param.use_l1
            error('Not implemented yet')
        end
        paramsflann.algorithm = 'kdtree';
        %TODO : optimize parameters in function of the number of
        %points
        paramsflann.checks = 32;
        paramsflann.trees = 1;
        % Use flann library
        [NN, DD] = flann_search(Xl', Xu', k, paramsflann);
        NN = transpose(NN);
        DD = transpose(DD);
    else
        %Built in matlab knn search
        if param.use_l1
            kdt = KDTreeSearcher(Xl, 'distance', 'cityblock');
            [NN, DD] = knnsearch(kdt, Xu, 'k', k ,  'Distance','cityblock');
        else
            kdt = KDTreeSearcher(Xl, 'distance', 'euclidean');
            [NN, DD] = knnsearch(kdt, Xu, 'k', k );
        end                   
        
    end


    % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]

    for ii = 1:N
        spi((ii-1)*k+1:ii*k) = repmat(ii, k, 1);
        spj((ii-1)*k+1:ii*k) = NN(ii, :);
        if param.weighted
            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = mean(DD(:)); end
                spv((ii-1)*k+1:ii*k) = exp(-DD(ii,:)/param.sigma);
            else
                if ~isfield(param, 'sigma'), param.sigma = mean(DD(:))^2; end
                spv((ii-1)*k+1:ii*k) = exp(-DD(ii,:).^2/param.sigma);
            end
        else
            spv((ii-1)*k+1:ii*k) = 1;
        end
        
    end
             

    %Actually create the sparse matrix from the 3-col values

    
    We = sparse(spj, spi, spv, M, N);

    W = [sparse(M,M) ,   We         ; ...
          We'        , sparse(N,N) ];
    
    %Fill in the graph structure
    G.W = W;
    G.coords = [Xl; Xu];
    G.We = We';
    G.de = sum(We)';
    if param.use_l1
        G.type = 'KNN classify l1';
    else
        G.type = 'KNN classify l2';
    end
    G = gsp_graph_default_parameters(G);
    
    
end

