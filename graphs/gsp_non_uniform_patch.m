function [ G ] = gsp_non_uniform_patch(N, n_patches, y_max, nn_params)
%GSP_NON_UNIFORM_PATCH Graph with nodes in patches with different density
%
%   Usage:  G = gsp_non_uniform_patch(N, n_patches, y_max, nn_params)
%           G = gsp_non_uniform_patch(N)
%           G = gsp_non_uniform_patch()
%
%   Input parameters:
%       N           : Number of nodes (default 200)
%       n_patches   : How many patches with interchanging density level?
%                     (default: 4)
%       y_max       : What is the width of the path? (default: 1)
%       nn_params   : Parameters for the nearest neighbours graph
%                     construction after the non-uniform sampling of the
%                     points. This structure is exactly the same as the one
%                     needed by gsp_nn_graph.
%
%   Output parameters:
%         G     : Graph structure.
%
%
%
%   Example:
%
%          G = gsp_non_uniform_patch(200, 5);
%          figure; gsp_plot_graph(G)
%          G = gsp_non_uniform_patch(100, 5, 0.2);
%          figure; gsp_plot_graph(G)
%
%see also: gsp_non_uniform
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_non_uniform_patch.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

% Author : Vassilis Kalofolias
 

if nargin < 1
    N = 200; 
end
if nargin < 2
    n_patches = 4;
end
if nargin < 3
    y_max = 1;
end
if nargin < 4
    nn_params = struct();
end

if ~isfield(nn_params, 'type'), nn_params.type = 'radius'; end
if ~isfield(nn_params, 'use_flann'), nn_params.use_flann = 0; end
if ~isfield(nn_params, 'center'), nn_params.center = 1; end
if ~isfield(nn_params, 'rescale'), nn_params.rescale = 1; end
if ~isfield(nn_params, 'k'), nn_params.k = 10; end
if ~isfield(nn_params, 'epsilon'), nn_params.epsilon = 0.01; end
if ~isfield(nn_params, 'use_l1'), nn_params.use_l1 = 0; end
if ~isfield(nn_params, 'target_degree'), nn_params.target_degree = 20; end;
if ~isfield(nn_params, 'symmetrize_type'), nn_params.symetrize_type = 'average'; end

    
    
x = zeros(N,1);

% Number of samples per patch
n_samples = ones(n_patches, 1);
% every two patches the sampling is 4 times more dense
n_samples(1:2:end) = 4;
% round to integers
n_samples = round(n_samples * N/sum(n_samples));
% make sure the sum is N by adding samples to last patch
n_samples(end) = n_samples(end) + N - sum(n_samples);
% indices of first samples of each patch
patch_first_sample = cumsum([1; n_samples]);

for p = 1:n_patches
    % find the x of the samples for this patch.
    x(patch_first_sample(p): patch_first_sample(p+1)-1) = (p-1) + rand(n_samples(p), 1);
end

x = lin_map(sort(x), [0, n_patches]);
% sample y uniformly from [0,1]
y = rand(size(x));
y = lin_map(y, [0, y_max]);

coords = [x, y];

G = gsp_nn_graph(coords, nn_params);
G.type = 'non_uniform_patch';
G.coords = coords;
end




