function [indx, indy, dist, Xo1, Xo2, epsilon, NN, Dist] = gsp_nn_distanz(X1, X2, param)
%GSP_NN_DISTANZ Compute the nearest neighboor distances
%   Usage : [indx, indy, dist] = gsp_nn_distanz( X1 );
%           [indx, indy, dist] = gsp_nn_distanz( X1, X2 );
%           [indx, indy, dist] = gsp_nn_distanz( X1, X2, param );
%           [indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(...)
%
%   Input parameters:
%       X1          : Input points 1
%       X2          : Input points 2
%       param       : Structure of optional parameters
%
%   Output parameters:
%       indx        : Indices over x
%       indy        : Indices over y
%       dist        : Distances
%       Xo1         : Points 1 after rescaling
%       Xo2         : Points 2 after rescaling
%       epsilon     : Radius of the ball (if the ball is used!)
%       NN          : Indices of closest neighbours of each node
%       Dist        : Sorted distances for each node
%
%   This function computes the nearest neighboors of Xin.
%
%   Additional parameters
%   ---------------------
%
%   * *param.type*      : ['knn', 'radius'] - the type of graph (default 'knn')
%   * *param.use_flann* : [0, 1] - use the FLANN library (default 0)
%   * *param.use_full*  : [0, 1] - Compute the full distance matrix and then
%     sparsify it (default 0)
%   * *param.flan_checks*: int - Number of checks for FLANN (default 256)
%     the higher the more precise, but the slower. Please consider the
%     following values: a) 32 not precise and fast, b) precise enought and
%     still fast, c) 4096 precise and may be slow.
%   * *param.nb_cores*  : int - Number of cores for FLANN (default 1)
%   * *param.center*    : [0, 1] - center the data
%   * *param.rescale*   : [0, 1] - rescale the data (in a 1-ball)
%   * *param.sigma*     : float  - the variance of the distance kernel
%   * *param.k*         : int - number of neighbors for knn
%   * *param.epsilon*   : float - the radius for the range search
%   * *param.use_l1*    : [0, 1] - use the l1 distance
%
%   See also: gsp_nn_graph
%

% Authors: Nathanael Perraudin, Vassilis Kalofolias, Johan Paratte

% Date: 8 September 2015, revised Nov 2015
% Testing: test_gsp_distanz

if nargin<3
    param = struct;
end

if nargin<2
    X2 = X1;
end

%Parameters
if ~isfield(param, 'type'), param.type = 'knn'; end
if ~isfield(param, 'use_flann'), param.use_flann = 0; end
if ~isfield(param, 'use_full'), param.use_full = 0; end
if ~isfield(param, 'center'), param.center = 1; end
if ~isfield(param, 'rescale'), param.rescale = 0; end
if ~isfield(param, 'k'), param.k = 10; end
if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
if ~isfield(param, 'use_cosine'), param.use_cosine = 0;
if ~isfield(param, 'target_degree'), param.target_degree = 0; end
if ~isfield(param, 'flann_nbcores'), param.flann_nbcores = 1; end
if ~isfield(param, 'flann_checks'), param.flann_checks = 256; end

Xo = transpose([X1, X2]);
[Nel, Nfeatures] = size(Xo);
n1 = size(X1, 2);
n2 = size(X2, 2);

% test if the binaries of flann are working
% TODO: better way to check if it works?
if param.use_flann
    try
        paramsflann.algorithm = 'kdtree';
        paramsflann.checks = 32;
        paramsflann.trees = 1;
        tmp = rand(100,10);
        [NN, Dist] = flann_search(tmp', tmp', 3, paramsflann);
    catch
        warning('Flann not available (either not compiled or Library not loaded) , going for the slow algorithm!')
        param.use_flann = 0;
    end
end


if strcmpi(param.type,'knn') && param.k>n1
    error('Not enough samples')
end



%Center the point cloud
if (param.center)
    Xo = Xo - repmat(mean(Xo), [Nel, 1]);
end

%Rescale the point cloud
if (param.rescale)
    bounding_radius = 0.5 * norm(max(Xo) - min(Xo));
    scale = nthroot(Nel, min(Nfeatures, 3)) / 10;
    Xo = Xo .* (scale / bounding_radius);
end


Xo1 = Xo(1:n1, :);
Xo2 = Xo((n1+1) : (n1+n2), :);

if ~(exist('KDTreeSearcher.m','file')==2) && ~param.use_flann
    param.use_full = 1;
    warning(['The Statistics and Machine Learning Toolbox is not availlable.',...
        'This function will not scale.'])
end

if param.use_cosine 
   warning(['Cosine distance nn graph uses an exhaustive search']);
end

if param.use_full
    
    D = gsp_distanz(transpose(Xo1),transpose(Xo2));
    switch param.type
        case 'knn'
            [Ds,ind] = sort(D,1,'ascend');
            indy = reshape(repmat(1:n2,param.k,1),[],1);
            indx = reshape(ind(1:param.k,:),[],1);
            dist = reshape(Ds(1:param.k,:),[],1);
            epsilon = nan;
        case 'radius'
            D(D>param.epsilon) = 0;
            [indx, indy, dist ] = find(D);
            epsilon = param.epsilon;
        otherwise
            error('Unknown type : allowed values are knn, radius');
    end
else
    
    switch param.type
        %Connect the k NN
        case 'knn'
            k = param.k;
            
            
            
            
            %Find kNN for each point in X (Using a kdtree)
            if param.use_flann
                if param.use_l1
                    flann_set_distance_type('manhattan');
                    paramsflann.checks = param.flann_checks;
                    paramsflann.cores = param.flann_nbcores;
                    paramsflann.algorithm = 'kdtree';
                    paramsflann.trees = 1;
                    [NN, Dist] = flann_search(transpose(Xo1), transpose(Xo2),...
                        k, paramsflann);
                    NN = transpose(NN);
                    Dist = transpose(Dist);
                else
                    flann_set_distance_type('euclidean');
                    %TODO : optimize parameters in function of the number of
                    %points
                    
                    %  paramsflann.target_precision = 0.9;
                    paramsflann.checks = param.flann_checks;
                    %  paramsflann.iterations = 5;
                    
                    paramsflann.cores = param.flann_nbcores;
                    % Use flann library
                    %             idx = flann_build_index(transpose(Xo1),struct('algorithm','kdtree','trees',1));
                    %             [NN, Dist] = flann_search(idx, transpose(Xo2),...
                    %                 k, paramsflann);
                    paramsflann.algorithm = 'kdtree';
                    paramsflann.trees = 1;
                    [NN, Dist] = flann_search(transpose(Xo1), transpose(Xo2),...
                        k, paramsflann);
                    NN = transpose(NN);
                    % Flann search return the distance squared. I do not know why
                    Dist = transpose(sqrt(Dist));
                end
            else
                %Built in matlab knn search
                if ~isreal(Xo)
                    Xo12 = [real(Xo1),imag(Xo1)];
                    Xo22 = [real(Xo2),imag(Xo2)];
                    if param.use_l1
                        kdt = KDTreeSearcher(Xo12, 'distance', 'cityblock');
                        [NN, Dist] = knnsearch(kdt, Xo22, 'k', k , ...
                            'Distance','cityblock');
                    else
                        kdt = KDTreeSearcher(Xo12, 'distance', 'euclidean');
                        [NN, Dist] = knnsearch(kdt, Xo22, 'k', k );
                    end
                else
                    if param.use_l1
                        kdt = KDTreeSearcher(Xo1, 'distance', 'cityblock');
                        [NN, Dist] = knnsearch(kdt, Xo2, 'k', k , ...
                            'Distance','cityblock');
                    else
                        kdt = KDTreeSearcher(Xo1, 'distance', 'euclidean');
                        [NN, Dist] = knnsearch(kdt, Xo2, 'k', k );
                    end
                end
            end
            
            
            % OLD CODE:
            % Create index matrices for the sparse W
            %         for ii = 1:n2
            %             indy((ii-1)*k+1:ii*k) = repmat(ii, k, 1);
            %             indx((ii-1)*k+1:ii*k) = transpose(NN(ii, :));
            %             dist((ii-1)*k+1:ii*k) = Dist(ii,1:end);
            %         end
            %
            
            % VAS: got rid of the for-loop.
            % Create index matrices for the sparse W
            indy = kron((1:n2)', ones(k, 1));
            indx = reshape(NN', k*n2, 1);
            dist = reshape(Dist', k*n2, 1);
            
            epsilon = nan;
            
            %Connect all the epsilon-closest NN
        case 'radius'
            %Create KDTree for fast NN computation
            
            if param.use_l1
                kdt = KDTreeSearcher(Xo1, 'distance', 'cityblock');
            else
                kdt = KDTreeSearcher(Xo1, 'distance', 'euclidean');
            end
            
            if param.use_flann
                warning('Flann is only used for k-NN search, not for radius search.');
            end
            
            if (param.target_degree == 0)
                epsilon = param.epsilon;
            else
                target_d = floor(param.target_degree);
                [NN, Dist] = knnsearch(kdt, Xo2, 'k', target_d);
                avg_d = mean(Dist);
                epsilon = avg_d(target_d);
            end
            %Find all neighbors at distance <= epsilon for each point in X
            if param.use_l1
                [NN, Dist] = rangesearch(kdt, Xo2, epsilon, ...
                    'Distance' ,'cityblock');
            else
                [NN, Dist] = rangesearch(kdt, Xo2, epsilon, ...
                    'distance', 'euclidean' );
            end
            
            % Create index matrices for the sparse W
            % VAS: got rid of the for-loop.
            indx = double(cell2mat(NN(:)')');
            dist = double(cell2mat(Dist(:)')');
            
            % number of NN for each node        EXAMPLE: [2 3 2]
            nNN = cellfun('length', NN);
            
            indy = zeros(length(indx) + 1, 1);
            % positions where indy changes value  EXAMPLE: [1 3 6 8]
            pos_y = cumsum([1; nNN]);
            % positions in indy that change value  EXAMPLE: [1 0 1 0 0 1 0 1]
            indy(pos_y(1:end)) = 1;
            % [1 0 1 0 0 1 0 1]  ->   [1 1 2 2 2 3 3 4]
            indy = cumsum(indy(1:end-1));
            
            
        otherwise
            error('Unknown type : allowed values are knn, radius');
    end
    
    
end

Xo1 = transpose(Xo1);
Xo2 = transpose(Xo2);

end