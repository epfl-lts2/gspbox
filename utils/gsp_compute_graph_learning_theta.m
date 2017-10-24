function [theta, theta_min, theta_max] = gsp_compute_graph_learning_theta(Z, k, geom_mean, is_sorted)

if nargin < 3 || isempty(geom_mean)
    geom_mean = 0;
end

if nargin < 4
    is_sorted = 0;
end

% Z is the zero-diagonal pairwise distance matrix between nodes
[theta_min, theta_max] = gsp_compute_theta_bounds(Z, geom_mean, is_sorted);
theta_min = theta_min(k);
theta_max = theta_max(k);

if k > 1
    theta = sqrt(theta_min * theta_max);
else
    theta = theta_min * 1.1;
end
