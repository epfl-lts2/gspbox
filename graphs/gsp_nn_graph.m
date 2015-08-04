function [ G ] = gsp_nn_graph( Xin, param )
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
%   Example:
%
%           P = gsp_pointcloud('bunny');
%           param.type = 'knn';
%           G = gsp_nn_graph(P, param);
%           gsp_plot_graph(G);
%
%   Additional parameters
%   ---------------------
%
%    param.type      : ['knn', 'radius']   the type of graph (default 'knn')
%    param.use_flann : [0, 1]              use the FLANN library
%    param.center    : [0, 1]              center the data
%    param.rescale   : [0, 1]              rescale the data (in a 1-ball)
%    param.sigma     : float               the variance of the distance kernel
%    param.k         : int                 number of neighbors for knn
%    param.epsilon   : float               the radius for the range search
%    param.use_l1    : [0, 1]              use the l1 distance
%    param.symetrize_type*: ['average','full'] symetrization type (default 'full')
%
%   See also: gsp_pointcloud
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_nn_graph.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

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
    if ~isfield(param, 'center'), param.center = 1; end
    if ~isfield(param, 'rescale'), param.rescale = 1; end
    if ~isfield(param, 'k'), param.k = 10; end
    if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
    if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
    if ~isfield(param, 'target_degree'), param.target_degree = 0; end;
    if ~isfield(param, 'symetrize_type'), param.symetrize_type = 'average'; end

    % test if the binaries of flann are working
    if param.use_flann
       try
            paramsflann.algorithm = 'kdtree';
            paramsflann.checks = 32;
            paramsflann.trees = 1;
            tmp = rand(100,10);
            [NN, D] = flann_search(tmp', tmp', 3, paramsflann);
       catch
            warning('Flann not compiled, going for the slow algorithm!')
            param.use_flann = 0;
       end
    end
    
    [N, d] = size(Xin);
   
    Xout = Xin;
    
    %Center the point cloud
    if (param.center)
        Xout = Xin - repmat(mean(Xin), [N, 1]);
    end
    
    %Rescale the point cloud
    if (param.rescale)
        bounding_radius = 0.5 * norm(max(Xout) - min(Xout));
        scale = nthroot(N, min(d, 3)) / 10;
        Xout = Xout .* (scale / bounding_radius);
    end
    
    switch param.type
        %Connect the k NN
        case 'knn'
            k = param.k;
     
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
                [NN, D] = flann_search(Xout', Xout', k+1, paramsflann);
                NN = transpose(NN);
                D = transpose(D);
            else
                %Built in matlab knn search
                if ~isreal(Xout)                   
                    Xout2 = [real(Xout),imag(Xout)];
                    if param.use_l1
                        kdt = KDTreeSearcher(Xout2, 'distance', 'cityblock');
                        [NN, D] = knnsearch(kdt, Xout2, 'k', k + 1, 'Distance','cityblock');
                    else
                        kdt = KDTreeSearcher(Xout2, 'distance', 'euclidean');
                        [NN, D] = knnsearch(kdt, Xout2, 'k', k + 1);
                    end                   
                else
                    if param.use_l1
                        kdt = KDTreeSearcher(Xout, 'distance', 'cityblock');
                        [NN, D] = knnsearch(kdt, Xout, 'k', k + 1,  'Distance','cityblock');
                    else
                        kdt = KDTreeSearcher(Xout, 'distance', 'euclidean');
                        [NN, D] = knnsearch(kdt, Xout, 'k', k + 1);
                    end                   
                end
            end
            

            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = mean(D(:)); end
            else
                if ~isfield(param, 'sigma'), param.sigma = mean(D(:))^2; end
            end
            
            % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
            for ii = 1:N
                spi((ii-1)*k+1:ii*k) = repmat(ii, k, 1);
                spj((ii-1)*k+1:ii*k) = NN(ii, 2:end)';
                if param.use_l1
                    spv((ii-1)*k+1:ii*k) = exp(-D(ii,2:end)/param.sigma);
                else
                    spv((ii-1)*k+1:ii*k) = exp(-D(ii,2:end).^2/param.sigma);
                end
            end
            

            
        %Connect all the epsilon-closest NN
        case 'radius'
            %Create KDTree for fast NN computation
            
            if param.use_l1
                kdt = KDTreeSearcher(Xout, 'distance', 'cityblock');
            else
                kdt = KDTreeSearcher(Xout, 'distance', 'euclidean');
            end
            
            if (param.target_degree == 0) 
                epsilon = param.epsilon;
            else
                target_d = floor(param.target_degree);
                [NN, D] = knnsearch(kdt, Xout, 'k', target_d);
                avg_d = mean(D);
                epsilon = avg_d(target_d);
            end
            tic;
            %Find all neighbors at distance <= epsilon for each point in X
            if param.use_l1
               [NN, D] = rangesearch(kdt, Xout, epsilon,'Distance' ,'cityblock');
            else
               [NN, D] = rangesearch(kdt, Xout, epsilon, 'distance', 'euclidean' );
            end
            toc;
            
            %Counting non-zero elements
            count = 0;
            for ii = 1:N
               count = count + length(NN{ii}) - 1; 
            end
            spi = zeros(count,1);
            spj = zeros(count,1);
            spv = zeros(count,1);
            start = 1;
            

            % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
            for ii = 1:N
                len = length(NN{ii}) - 1;
                spi(start:start+len-1) = repmat(ii, len, 1);
                spj(start:start+len-1) = NN{ii}(2:end)';
                if param.use_l1
                    if ~isfield(param, 'sigma'), param.sigma = param.epsilon/5; end
                    spv(start:start+len-1) = exp(-D{ii}(2:end)/param.sigma);
                else
                    if ~isfield(param, 'sigma'), param.sigma = param.epsilon/5; end
                    spv(start:start+len-1) = exp(-D{ii}(2:end).^2/param.sigma);
                end
                start = start + len;
            end
            
        otherwise
            error('Unknown type : allowed values are knn, radius');
    end

    %Actually create the sparse matrix from the 3-col values
    W = sparse(spi, spj, spv, size(Xin,1), size(Xin,1));
    
    % Computes the average degree when using the epsilon-based neighborhood
    if (strcmp(param.type,'radius'))
        text = sprintf('Average degree = %d', nnz(W)/size(W, 1));
        disp(text);
    end

    % Sanity check
    if size(W,1) ~= size(W,2), error('Weight matrix W is not square'); end
    
    % Symmetry checks
    %if issymmetric(W)
    if (norm(W - W', 'fro') == 0)
        disp('The matrix W is symmetric');
    else
         W = gsp_symetrize(W,param.symetrize_type);
    end
    
    %Fill in the graph structure
    G.N = N;
    G.W = W;
    G.coords = Xout;
    %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
    if param.use_l1
        G.type = 'nearest neighbors l1';
    else
        G.type = 'nearest neighbors';
    end
    %G.vertex_size=30;

    G = gsp_graph_default_parameters(G);
end


