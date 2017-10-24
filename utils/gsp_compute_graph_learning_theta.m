function [theta, theta_min, theta_max] = gsp_compute_graph_learning_theta(Z, k, geom_mean, is_sorted)

if nargin < 3 || isempty(geom_mean)
    geom_mean = 0;
end

if nargin < 4
    is_sorted = 0;
end

% Z is the zero-diagonal pairwise distance matrix between nodes
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/utils/gsp_compute_graph_learning_theta.html

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
[theta_min, theta_max] = gsp_compute_theta_bounds(Z, geom_mean, is_sorted);
theta_min = theta_min(k);
theta_max = theta_max(k);

if k > 1
    theta = sqrt(theta_min * theta_max);
else
    theta = theta_min * 1.1;
end

