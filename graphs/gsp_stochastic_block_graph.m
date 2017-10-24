function G = gsp_stochastic_block_graph(N, k, params)
%GSP_STOCHASTIC_BLOCK_GRAPH  Create a stochastic block graph
%   Usage:  G = gsp_stochastic_block_graph( N );
%           G = gsp_stochastic_block_graph(N , k);
%
%   Input parameters:
%         N     : Number of nodes (default 1024)
%         k     : Number of clusters (default 5)
%         params: Structure of optional parameters
%   Output parameters:
%         G     : Graph structure.
%
%   param is an optional structure with the following fields
%
%    params.p        : Intra-cluster edge probability (default 0.7)
%    params.q        : Inter-cluster edge probability (default 0.3/k)
%    params.z        : Assignment vector of nodes (default uniform random)
%    params.M        : Link probability matrix between clusters (default uses p and q)
%    params.directed : Flag the graph as directed or not (default false)
%
%   Use the stochastic block model to create a graph.
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_stochastic_block_graph.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

% Author: Pierre Vandergheynst, Nathanael Perraudin
% Date  : 2 novemeber 2015 (revision: 6 october 2016 -- Lionel Martin)



%% Stochastic Block Model generator
if nargin<1
    N = 1024; % number of nodes
end
if nargin<2
    k = 5; % number of clusters
end
if nargin < 3
    params = struct;
end

if ~isfield(params, 'p')
    params.p = 0.7;
end

if ~isfield(params, 'q')
    params.q = (1-params.p) / k;
end

if ~isfield(params, 'force_full')
    params.force_full = 0;
end

if ~isfield(params, 'auto_gen_M')
    params.auto_gen_M = 0;
end

if (~isfield(params, 'z') || length(params.z) ~= N || max(params.z) > k || min(params.z) < 1)
    params.z = randi(k, 1, N);
end

if (~isfield(params, 'M') || size(params.M) ~= [k, k])
    if params.force_full
        params.M = params.q * ones(k);
        params.M(1:k+1:end) = params.p;
    end
    params.auto_gen_M = 1;
end

if ~isfield(params, 'directed')
    params.directed = false;
end

%% partition with clusters of homogenous sizes
%
% THIS PART IS NOT NECESSARY BECAUSE SLOW AND REPLACED
% BY THE SIMPLE "UNIFORM" CLUSTER ASSIGNMENT ABOVE
%
% L = ceil(N/k);
% for i=1:k-1
%     z((i-1)*L + 1:i*L) = i;
% end
% z((k-1)*L + 1:end) = k;
 
%% Generate adjacency
z = params.z;

% for i=1:N
%     for j=i+1:N
%         W(i, j) = ( rand <= M(z(i), z(j)) );
%         if params.directed
%             W(j, i) = ( rand <= M(z(j), z(i)) );
%         else
%             W(j, i) = W(i, j);
%         end
%     end
% end
if params.auto_gen_M && ~params.force_full
    [val_z, idx_z] = sort(z);
    counts = diff(find(diff([0, val_z, N])));

    if length(counts) ~= k
        error('There is at least one empty class. Check your z.');
    end

    W = sprandsym(counts(k), params.p);
    nb_cols_rect = 0;

    for i=k:-1:2
        nb_cols_rect = nb_cols_rect + counts(i);
        rect = sprand(counts(i-1), nb_cols_rect, params.q);
        top_left = sprandsym(counts(i-1), params.p);
        W = vertcat(horzcat(top_left, rect), horzcat(rect', W));
    end

    W(1:N+1:end) = 0;
    G.W(idx_z, idx_z) = abs(W) > 0;

else
    M = params.M;

    if params.directed
        W = rand(N) <= M(z, z);
    else
        A = rand(N);
        A(logical(triu(ones(N)))) = 1;
        W = A <= M(z, z);
        W = W + W';
    end

    G.W = sparse(W);
end
 
G = gsp_graph_default_parameters(G);
G.info.node_com = z;

G.coords = ones(N, 2);
com_coords = sqrt(N) * [-cos(2*pi*(1:k)/k)', sin(2*pi*(1:k)/k)'];

% create uniformly random points in the unit disc
for ii = 1:N
    % use rejection sampling to sample from a unit disc (probability = pi/4)
    while norm(G.coords(ii, :)) >= 1/2
        % sample from the square and reject anything outside the circle
        G.coords(ii, :) = [rand-.5, rand-.5];
    end
end

% add the offset for each node depending on which community it belongs to
for ii = 1:k
    idx_ii = find(z==ii);
    rad_com = sqrt(numel(idx_ii));
    G.coords(idx_ii, :) = bsxfun(@plus, rad_com * G.coords(idx_ii, :), com_coords(ii, :));
end

end

