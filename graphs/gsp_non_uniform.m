function [ G ] = gsp_non_uniform(N, r)
%GSP_NON_UNIFORM Create a random graph from non_uniform sampling of points
%   Usage:  G = gsp_non_uniform(N, r)
%           G = gsp_non_uniform(N)
%           G = gsp_non_uniform()
%
%   Input parameters:
%         N     : Number of nodes (default 200)
%         r     : Factor of non-uniformity (default 3)
%                 the histogram of x axis of samples is proportional to
%                 exp(x*r)
%   Output parameters:
%         G     : Graph structure.
%
%
%   Example:::
%
%          G = gsp_non_uniform(200, 1);
%          figure; gsp_plot_graph(G)
%          G = gsp_non_uniform(200, 5);
%          figure; gsp_plot_graph(G)
%
%
% see also: gsp_sensor, gsp_non_uniform_patch

% Author : Vassilis Kalofolias
 

% Optional input arguments
% if nargin < 2 
%     k = 6;
% end

if nargin < 1
    N = 200; 
end
if nargin < 2
    r = 3;
end

% sample x exponentially from [0,5]
x = exp(rand(N, 1) * r);
x = lin_map(sort(x), [0, 5]);
% sample y uniformly from [0,1]
y = rand(N, 1);
y = lin_map(y, [0, 1]);

coords = [x, y];

D = gsp_distanz(coords');

Ds = sort(D, 'ascend');
%max_min_dist = max(min(D + realmax * eye(N)));
s2 = mean(Ds(2,:).^2);
%s2 = mean(Ds(5,:).^2) / 9;

W = exp(- D.^2 / s2 / 4);
W(1:N+1:end) = 0;

% keep a neighbour for the most distant 
min_kept_w = min(max(W));

% 


% Kill all small connections but always have at least one neighbor
W = W .* (W >= min(min_kept_w, 0.1));

% binary?
%W = double(W>0);

G = gsp_graph(W, coords, [0 5 0 1]);
G.type = 'non_uniform';

end



