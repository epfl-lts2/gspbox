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
%   Example:
%
%           P = rand(100,2);
%           G = gsp_nn_hypergraph(P)
%
%   Additional parameters
%   ---------------------
%
%    param.use_flann : [0, 1]              use the FLANN library
%    param.center    : [0, 1]              center the data
%    param.rescale   : [0, 1]              rescale the data (in a 1-ball)
%    param.sigma     : float               the variance of the distance kernel
%    param.k         : int                 number of neighbors for knn
%
%   See also: gsp_nn_graph
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_nn_hypergraph.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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


