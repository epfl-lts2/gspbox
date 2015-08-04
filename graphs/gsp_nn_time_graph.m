function [ G ] = gsp_nn_time_graph( Xin, Nf, param )
%GSP_PCL_NN_GRAPH Create a nearest neighbors graph from a point cloud
%   Usage :  G = gsp_nn_time_graph( Xin, Nf );
%            G = gsp_nn_time_graph( Xin, Nf, param );
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
%    param.symetrize_type*: ['average','full'] symetrization type (default 'full')
%
%   See also: gsp_pointcloud
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_nn_time_graph.php

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

% Author: Johan Paratte
% Date: 16 June 2014
% 

    if nargin < 2
    % Define parameters
        param = {};
    end
    
    %Parameters
    if ~isfield(param, 'type'), param.type = 'knn'; end
    if ~isfield(param, 'use_flann'), param.use_flann = 0; end
    if ~isfield(param, 'center'), param.center = 1; end
    if ~isfield(param, 'rescale'), param.rescale = 1; end
    if ~isfield(param, 'k_in'), param.k_in = 10; end
    if ~isfield(param, 'k_out'), param.k_out = 5; end
    if ~isfield(param, 'sigma'), param.sigma = 0.1; end
    if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
    if ~isfield(param, 'symetrize_type'), param.symetrize_type = 'average'; end

    
    %[M, d] = size(Xin);
    
%     if ~(mod(M, Nf) == 0)
%         error('Number of frames not coherent with input size');
%     end
    
    %N = M / Nf;
    
    Ntotal = 0;
    dim = 0;
    for n = 1:Nf
        [N, d] = size(Xin{n});
        dim = d;
        Ntotal = Ntotal + N;
    end
    
    Xout = Xin;
    
    %Center the point cloud
    if (param.center)
        for n = 1:Nf
            [N, d] = size(Xin{n});
            Xout{n} = Xin{n} - repmat(mean(Xin{n}), [N, 1]);
        end
    end
    
    %Rescale the point cloud
    if (param.rescale)
        for n = 1:Nf
            [N, d] = size(Xin{n});
            bounding_radius = 0.5 * norm(max(Xout{n}) - min(Xout{n}));
            scale = nthroot(N, min(d, 3)) / 10;
            Xout{n} = Xout{n} .* (scale / bounding_radius);
        end
    end
    
    switch param.type
        %Connect the k NN
        case 'knn'
            ki = param.k_in;
            ko = param.k_out;
            
            nnz = 0;
            
            for n = 1:Nf
               [N, d] = size(Xin{n});
               nnz = nnz + N*ki;
            end
            
            for n = 1:Nf-1
               [N, d] = size(Xin{n});
               nnz = nnz + N*ko;
            end
            
            nnz
            
            spi = zeros(nnz,1);
            spj = zeros(nnz,1);
            spv = zeros(nnz,1);
            
            offset = 0;
            linear_offset = 0;
            
            for n = 1:Nf
                [N, d] = size(Xin{n});
                %Find kNN for each point in X (Using a kdtree)
                if param.use_flann
                    paramsflann.algorithm = 'kdtree';
                    %TODO : optimize parameters in function of the number of
                    %points
                    paramsflann.checks = 32;
                    paramsflann.trees = 1;
                    % Use flann library
                    [NN, D] = flann_search(Xout{n}', Xout{n}', ki+1, paramsflann);
                    NN = transpose(NN);
                    D = transpose(D);
                else
                    %Built in matlab knn search
                    if ~isreal(Xout{n})                   
                        Xout2 = [real(Xout{n}),imag(Xout{n})];
                        kdt = KDTreeSearcher(Xout2, 'distance', 'euclidean');
                        [NN, D] = knnsearch(kdt, Xout2, 'k', ki + 1);
                    else
                        kdt = KDTreeSearcher(Xout{n}, 'distance', 'euclidean');
                        [NN, D] = knnsearch(kdt, Xout{n}, 'k', ki + 1);
                    end
                end

                % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
                for ii = 1:N
                    spi(offset+((ii-1)*ki+1):offset+(ii*ki)) = repmat(linear_offset + ii, ki, 1);
                    spj(offset+((ii-1)*ki+1):offset+(ii*ki)) = NN(ii, 2:end)' + repmat(linear_offset, ki, 1);
                    spv(offset+((ii-1)*ki+1):offset+(ii*ki)) = exp(-D(ii,2:end).^2/param.sigma);
                end
                
                clear NN
                clear D
                clear kdt
                
                offset = offset + N*ki
                linear_offset = linear_offset + N
            end
            
            linear_offset = 0;
            for n = 1:Nf-1
                [N, d] = size(Xin{n});
                %Find kNN for each point in X (Using a kdtree)
                if param.use_flann
                    paramsflann.algorithm = 'kdtree';
                    %TODO : optimize parameters in function of the number of
                    %points
                    paramsflann.checks = 32;
                    paramsflann.trees = 1;
                    % Use flann library
                    [NN, D] = flann_search(Xout{n+1}', Xout{n}', ko+1, paramsflann);
                    NN = transpose(NN);
                    D = transpose(D);
                else
                    %Built in matlab knn search
                    if ~isreal(Xout{n})                   
                        Xout2 = [real(Xout{n+1}),imag(Xout{n+1})];
                        kdt = KDTreeSearcher(Xout2, 'distance', 'euclidean');
                        Xout2 = [real(Xout{n}),imag(Xout{n})];
                        [NN, D] = knnsearch(kdt, Xout2, 'k', ko + 1);
                    else
                        kdt = KDTreeSearcher(Xout{n+1}, 'distance', 'euclidean');
                        [NN, D] = knnsearch(kdt, Xout{n}, 'k', ko + 1);
                    end
                end

                % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
                for ii = 1:N
                    spi(offset+((ii-1)*ko+1):offset+(ii*ko)) = repmat(linear_offset + ii, ko, 1);
                    spj(offset+((ii-1)*ko+1):offset+(ii*ko)) = NN(ii, 2:end)' + repmat(linear_offset+N, ko, 1);
                    spv(offset+((ii-1)*ko+1):offset+(ii*ko)) = exp(-D(ii,2:end).^2/param.sigma);
                end
                
                clear NN
                clear D
                clear kdt
                
                offset = offset + N*ko
                linear_offset = linear_offset+N
            end
            
            %Actually create the sparse matrix from the 3-col values
            W = sparse(spi, spj, spv, Ntotal, Ntotal);
            
        %Connect all the epsilon-closest NN
        case 'radius'
            error('Radius mode non supported yet');
            %Create KDTree for fast NN computation
            kdt = KDTreeSearcher(Xout, 'distance', 'euclidean');
            
            eps = param.epsilon;
            tic;
            %Find all neighbors at distance <= epsilon for each point in X
            [NN, D] = rangesearch(kdt, Xout, eps);
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
            
            tic;
            % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
            for ii = 1:N
                len = length(NN{ii}) - 1;
                spi(start:start+len-1) = repmat(ii, len, 1);
                spj(start:start+len-1) = NN{ii}(2:end)';
                spv(start:start+len-1) = exp(-D{ii}(2:end).^2/param.sigma);
                start = start + len;
            end
            toc;
            tic;
            %Actually create the sparse matrix from the 3-col values
            W = sparse(spi, spj, spv, size(Xin,1), size(Xin,1));
        otherwise
            error('Unknown type : allowed values are knn, radius');
    end
    
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
    coords = zeros(Ntotal, dim);
    offset = 0;
    for n =1:Nf
        [N, d] = size(Xin{n});
        coords(offset + 1:offset + N, :) = Xin{n};
        offset = offset + N;
    end
    G.coords = coords;
    %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
    G.type = 'nearest neighbors';
    %G.vertex_size=30;

    G = gsp_graph_default_parameters(G);
end


