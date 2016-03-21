function [f_interpolated]=gsp_interpolate(G,f_subsampled,keep_inds,param)
%GSP_INTERPOLATE Interpolation of a graph signal
%   Usage:  f_interpolated=gsp_interpolate(f_subsampled, G, keep_inds);
%           f_interpolated=gsp_interpolate(f_subsampled, G, keep_inds, param);
%
%   Input parameters:
%         f_subsampled     : A signal on the subset of the vertices of G indexed by keep_inds
%         G                : Graph structure.
%         keep_inds        : The vertex indices of V_1
%   Output parameters:
%         f_interpolated   : Interpolated graph signal on G.
%
%   'gsp_interpolate(f_subsampled,G,keep_inds)' performs the interpolation:
%
%      f_interpolated = \Phi_{V_1} \Phi_{V_1,V_1}^{-1} f_{V_1} 
%                     = \Phi_{V_1} [\bar{L}_{V_1,V_1}-\bar{L}_{V_1,V_1^c}(\bar{L}_{V_1^c,V_1^c})^{-1}\bar{L}_{V_1^c,V_1}] f_{V_1} 
%
%   Note that keep_inds are the vertex indices of V_1; i.e., those
%   where the signal values are available to use in the interpolation. 
%
%   param contains the following fields
%
%    param.order : Degree of the Chebyshev approximation (default=100 
%     to yield a good approximation of the default Green's kernel). This
%     parameter is passed to gsp_filter_analysis.
%    param.regularize_epsilon  : The regularized graph Laplacian is
%     bar{L}=L+epsilon I (default = 0.005). A smaller epsilon may 
%     lead to better regularization, but will also require a higher order 
%     Chebyshev approximation. 
%
%   See also: gsp_filter_analysis 
% 
%   References:
%     I. Pesenson. Variational splines and paley--wiener spaces on
%     combinatorial graphs. Constructive Approximation, 29(1):1--21, 2009.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_interpolate.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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

%   See also: gsp_resgression_tik

%   Author : David I Shuman, Nathanael Perraudin
%   Date   : 26 November 2015
%   Testing: test_pyramid
  
% Read input parameters
if nargin < 4
    param = struct;
end
    
if ~isfield(param,'order'), param.order = 100; end 
if ~isfield(param, 'regularize_epsilon'), param.regularize_epsilon=.005; end

% Compute coefficients alpha for each of the Green's functions translated
% to center vertices in V_1
elim_inds=setdiff(1:G.N,keep_inds);
regularized_L= G.L+param.regularize_epsilon*eye(G.N);
alpha_upsampled=zeros(G.N,size(f_subsampled,2));
alpha_upsampled(keep_inds,:) = ...
    regularized_L(keep_inds,keep_inds) * f_subsampled ...
    -regularized_L(keep_inds,elim_inds) * ...
    ( regularized_L(elim_inds,elim_inds) \ ...
    ( regularized_L(elim_inds,keep_inds) * f_subsampled));
    
% Now compute f_interp=\Phi_{V_1} \alpha
green_kernel=@(x) 1./(x+param.regularize_epsilon);
f_interpolated=gsp_filter_analysis(G, green_kernel ,alpha_upsampled,param);

end
    
