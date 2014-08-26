function [ G ] = gsp_pcl_nn_graph( Xin, param )
%GSP_PCL_NN_GRAPH Create a nearest neighbors graph from a point cloud
%   Usage :  G = gsp_pcl_nn_graph( Xin );
%            G = gsp_pcl_nn_graph( Xin, param );
%
%   Input parameters:
%       Xin         : Input points
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_pcl_nn_graph( Xin, param )' creates a graph from positional data. The points are 
%   connected to their neighbors (either belonging to the k nearest 
%   neighbors or to the epsilon-closest neighbors. 
%
%   Example:
%
%           P = gsp_pointcloud('bunny');
%           param.type = 'knn';
%           G = gsp_pcl_nn_graph(P, param);
%           gsp_plot_graph(G);
%
%   Additional parameters
%   ---------------------
%
%    param.type      : ['knn', 'radius']   the type of graph
%    param.center    : [0, 1]              center the data
%    param.rescale   : [0, 1]              rescale the data (in a 1-ball)
%    param.sigma     : float               the variance of the distance kernel
%    param.k         : int                 number of neighbors for knn
%    param.epsilon   : float               the radius for the range search
%
%   See also: gsp_pointcloud
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_pcl_nn_graph.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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
    if ~isfield(param, 'center'), param.center = 1; end
    if ~isfield(param, 'rescale'), param.rescale = 1; end
    if ~isfield(param, 'k'), param.k = 10; end
    if ~isfield(param, 'sigma'), param.sigma = 0.1; end
    if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end

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
    
    %Create KDTree for fast NN computation
    kdt = KDTreeSearcher(Xout, 'distance', 'euclidean');
    
    switch param.type
        %Connect the k NN
        case 'knn'
            k = param.k;
            spidx = zeros(N*k, 3);
            
            %Find kNN for each point in X
            [NN, D] = knnsearch(kdt, Xout, 'k', k + 1);
            
            % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
            for ii = 1:N
                spidx((ii-1)*10+1:ii*10, 1) = repmat(ii, k, 1);
                spidx((ii-1)*10+1:ii*10, 2) = NN(ii, 2:end)';
                spidx((ii-1)*10+1:ii*10, 3) = exp(-D(ii,2:end).^2/param.sigma);
            end
            
            %Actually create the sparse matrix from the 3-col values
            W = spconvert(spidx);
            
        %Connect all the epsilon-closest NN
        case 'radius'
            eps = param.epsilon;
            %Find all neighbors at distance <= epsilon for each point in X
            [NN, D] = rangesearch(kdt, Xout, eps);
            
            %Counting non-zero elements
            count = 0;
            for ii = 1:N
               count = count + length(NN{ii}) - 1; 
            end
            spidx = zeros(count, 3);
            start = 1;
            
            % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
            for ii = 1:N
                len = length(NN{ii}) - 1;
                spidx(start:start+len-1, 1) = repmat(ii, len, 1);
                spidx(start:start+len-1, 2) = NN{ii}(2:end)';
                spidx(start:start+len-1, 3) = exp(-D{ii}(2:end).^2/param.sigma);
                start = start + len;
            end
            
            %Actually create the sparse matrix from the 3-col values
            W = spconvert(spidx);
            
        otherwise
            error('Unknown type : allowed values are knn, radius');
    end
    
    % Computes the average degree when using the epsilon-based neighborhood
    if (strcmp(param.type,'radius'))
        text = sprintf('Average degree = %d', nnz(W)/size(W, 1));
        disp(text);
    end

    % Symmetry checks
    if (norm(W - W', 'fro') == 0)
        disp('The matrix W is symmetric');
    else
        disp(['WARNING : The matrix W is not symmetric !',...
              ' Artificial symmetrization']);
        W = (W+W')/2;
    end
    
    %Fill in the graph structure
    G.N = N;
    G.W = W;
    G.coords = Xout;
    %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
    G.type = 'nearest neighbors';
    %G.vertex_size=30;

    G = gsp_graph_default_parameters(G);
end


