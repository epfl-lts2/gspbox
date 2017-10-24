% Compute the values of parameter theta (controlling sparsity) that should
% be expected to give each sparsity level. Return upper and lower bounds 
% each sparsity level k=[1, ..., n-1] neighbors/node.
%
% USAGE: [theta_l, theta_u, Z_sorted] = gsp_compute_theta_bounds(Z)
%                                       gsp_compute_theta_bounds(Z, geom_mean)
%                                       gsp_compute_theta_bounds(Z, geom_mean, is_sorted)
% 
% geom_mean:    use geometric mean instead of arithmetic mean? default: 0
% is_sorted:    is Z already sorted? default: 0
%
%
%
% 
%
% code author: Vassilis Kalofolias
% date: Aug 2016
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/utils/gsp_compute_theta_bounds.html

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


function [theta_l, theta_u, Z_sorted] = gsp_compute_theta_bounds(Z, geom_mean, is_sorted)

if nargin < 2 || isempty(geom_mean)
    geom_mean = 0;
end

if nargin < 3
    is_sorted = 0;
end

% Z is the zero-diagonal pairwise distance matrix between nodes


if ismatrix(Z)
    
    if is_sorted
        Z_sorted = Z;
    else
        n = length(Z);
        % don't take into account the diagonal of Z
        Z_sorted = zeros(n, n-1);
        for i=1:n
            Z_sorted(i, :) = sort(Z(i, [1:i-1, i+1:n]), 2, 'ascend');
        end
    end
    [m, n] = size(Z_sorted);
    
    B_k = cumsum(Z_sorted, 2);       % cummulative sum for each row
    K_mat = repmat((1:n), m, 1);

    %% Theoretical intervals of theta for each desired sparsity level:
    if geom_mean == 0
        theta_u = mean(1./sqrt(K_mat.*Z_sorted.^2 - B_k.*Z_sorted));
    else
    % try geometric mean instead of arithmetic:
        theta_u = exp(mean(log(1./sqrt(K_mat.*Z_sorted.^2 - B_k.*Z_sorted))));
    end 
    theta_l = [theta_u(2:end), 0];
end




