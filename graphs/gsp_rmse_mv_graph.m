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
%   nan value. But it is much slower than GSP_NN_GRAPH.
%
%   Additional parameters
%   ---------------------
%
%    param.type      : ['knn', 'radius']   the type of graph (default 'knn')
%    param.sigma     : float               the variance of the distance kernel
%    param.k         : int                 number of neighbors for knn
%    param.epsilon   : float               the radius for the range search
%    param.symetrize_type*: ['average','full'] symetrization type (default 'full')
%    param.center    : [0, 1]              center the data
%    param.rescale   : [0, 1]              rescale the data (in a 1-ball)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_rmse_mv_graph.php

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

% Author: Nathanael Perraudin
% Date : 12 March 2015
% Testing: test_rmse

if nargin<2
    param = struct;
end

if ~isfield(param, 'sigma'), param.sigma = 1; end
if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'k'), param.k = 10; end
if ~isfield(param, 'symetrize_type'), param.symetrize_type = 'average'; end
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
C = gsp_rmse_mv(Xout',p);
C = C*sqrt(d);
sigma = param.sigma;


% sparsification
switch param.type
    case 'knn'
        W = exp(-C.^2/sigma);
        W = W-diag(diag(W));
        for ii=1:size(W,1)
            [~,sortIndex] = sort(W(ii,:),'descend');  %# Sort the values in
                                                      %#   descending order
            W(ii,sortIndex(param.k+1:end)) = 0;
        end
        W = gsp_symetrize(W,param.symetrize_type);
    case 'radius'
        C(C>param.epsilon) = 0;
        W = exp(-C.^2/sigma);
        W = W-diag(diag(W));
    otherwise
       error('Unknown type : allowed values are knn, radius');
end


% Create the graph
G = gsp_graph(sparse(W));


        
end

