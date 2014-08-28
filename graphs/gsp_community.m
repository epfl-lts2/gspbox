function [G] = gsp_community(N, param)
%GSP_COMMUNITY Create a community graph
%   Usage: G = gsp_community(N);
%          G = gsp_community();
%          G = gsp_community(N, param );
%
%   Input parameters
%       - N     : Number of nodes (default 256)
%       - param : Structure of optional parameters
%   Output parameters
%       - G     : Graph
%
%   This function create a 2 dimentional random sensor graph. All the
%   coordonates are between 0 and 1.
%   
%   param is an optional structure with the following fields
%
%    param.Nc : Number of communities (default round(sqrt(N)/2) )
%    param.verbose*: display parameter - 0 no log - 1 display the errors
%     (default 1) 
%    param.com_sizes : size of the communities. The sum of the sizes has
%     to be equal to N. Leave this field empty if you want random sizes.
%    param.min_comm : Minimum size of the community 
%     (default round(N / param.Nc / 3) )
%    param.min_deg: Minimum degree of each nodes 
%     (default round(param.min_comm/2)) (NOT WORKING YET!)
%    param.size_ratio*: ratio between radius of world and radius of
%     communities (default 1)
%    param.world_density  probability of a random edge between any pair
%     of nodes (default 1/N)
%
%   Example:
%
%          G = gsp_community();
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_community.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

% Author: Vassilis Kalofolias, Nathanael Perraudin
% Date: March 2014
% Testing: test_graphs


% TODO: This function has to be revisited

if nargin < 2
   param = struct;
end
if nargin < 1
   N = 256; 
end

if ~isfield(param, 'Nc'), param.Nc = round(sqrt(N)/2); end
if isfield(param, 'com_sizes')
    if sum(param.com_sizes) ~= N
        error(['GSP_COMMUNITY: The sum of the community sizes has ',...
            'to be equal to N']);
    end
else
    param.com_sizes = []; 
end
if ~isfield(param, 'min_comm'), param.min_comm = round(N / param.Nc / 3); end
if ~isfield(param, 'min_deg'), param.min_deg = round(param.min_comm/2); end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'param.size_ratio'), param.size_ratio = 1; end     
if ~isfield(param, 'world_density'), param.world_density = 1/N; end     



if isempty(param.com_sizes)
    com_lims = sort(randperm(N - (param.min_comm-1) * param.Nc - 1, param.Nc-1), 'ascend');
    com_lims = com_lims + cumsum((param.min_comm-1) * ones(size(com_lims)));
    com_lims = [0, com_lims, N];
    param.com_sizes = diff(com_lims);
end

if param.verbose > 2
    X = zeros(10000, param.Nc + 1);
    %pick randomly param.Nc-1 points to cut the rows in communtities:
    for i=1:10000
        com_lims_temp = sort(randperm(N - (param.min_comm-1) * param.Nc - 1, param.Nc-1), 'ascend');
        com_lims_temp = com_lims_temp + cumsum((param.min_comm-1) * ones(size(com_lims_temp)));
        X(i,:) = [0, com_lims_temp, N];
    end
    dX = diff(X')';
    for i=1:param.Nc; figure;hist(dX(:,i), 100); title('histogram of row community size'); end
    clear X com_lims_temp
end


rad_world = param.size_ratio * sqrt(N);
com_coords = rad_world * [-cos(2*pi*(1:param.Nc)/param.Nc)', sin(2*pi*(1:param.Nc)/param.Nc)'];

G.coords = ones(N, 2);

% create uniformly random points in the unit disc
for i = 1:N
    % use rejection sampling to sample from a unit disc (probability = pi/4)
    while norm(G.coords(i, :)) >= 1/2
        % sample from the square and reject anything outside the circle
        G.coords(i, :) = [rand-.5, rand-.5];
    end
end

% add the offset for each node depending on which community it belongs to
info.node_com = zeros(N, 1);
for i = 1:param.Nc
    com_size = param.com_sizes(i);
    rad_com = sqrt(com_size);

    node_ind = (com_lims(i) + 1) : com_lims(i+1);
    G.coords(node_ind, :) = bsxfun(@plus, rad_com * G.coords(node_ind, :), com_coords(i, :));
    info.node_com(node_ind) = i;
end
    

% TODO: this can (and should to prevent overlap) be done for each community separately!
D = gsp_distanz(G.coords');
%R = rmse_mv(G.coords')*sqrt(2);
%W = graph_k_NN(exp(-D.^2), param.min_deg);
W = exp(-D.^2);
W(W<1e-3) = 0;


%TODO: this could be more sophisticated (e.g. one sigma for communities,
%one sigma for inter-community connections, exp(-d^2/sigma) weights!
W = W + abs(sprandsym(N, param.world_density));
W = double(abs(W) > 0);

G.W = sparse(W);
G.type = 'Community';

% return additional info about the communities
info.com_lims = com_lims;
info.com_coords = com_coords;
info.com_sizes = param.com_sizes;

G.info = info;
G = gsp_graph_default_parameters( G );


end






