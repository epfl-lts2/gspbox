function [f_interpolated]=gsp_interpolate(G,f_subsampled,keep_inds,param)
%GSP_INTERPOLATE Interpolation of a graph signal
%   Usage:  f_interpolated=gsp_interpolate(f_subsampled, G, keep_inds);
%           f_interpolated=gsp_interpolate(f_subsampled, G, keep_inds, param);
%
%   Input parameters:
%         f_subsampled     : A signal on the subset of the vertices of G indexed by keep_inds
%         G                : Graph structure.
%         keep_inds        : The vertex indices of $V_1$
%   Output parameters:
%         f_interpolated   : Interpolated graph signal on G.
%
%   'gsp_interpolate(f_subsampled,G,keep_inds)' performs the interpolation:
%
%   .. f_interpolated = \Phi_{V_1} \Phi_{V_1,V_1}^{-1} f_{V_1} 
%   ..                = \Phi_{V_1} [\bar{L}_{V_1,V_1}-\bar{L}_{V_1,V_1^c}(\bar{L}_{V_1^c,V_1^c})^{-1}\bar{L}_{V_1^c,V_1}] f_{V_1} 
%
%   .. math:: f_{interpolated}=\Phi_{{\cal V}_1} \Phi_{{\cal V}_1,{\cal V}_1}^{-1} f_{{\cal V}_1} = \Phi_{{\cal V}_1} \left[\bar{{\cal L}}_{{\cal V}_1,{\cal V}_1}-\bar{{\cal L}}_{{\cal V}_1,{\cal V}_1^c}(\bar{{\cal L}}_{{\cal V}_1^c,{\cal V}_1^c})^{-1}\bar{{\cal L}}_{{\cal V}_1^c,{\cal V}_1}\right] f_{{\cal V}_1}
%
%   Note that *keep_inds* are the vertex indices of $V_1$; i.e., those
%   where the signal values are available to use in the interpolation. 
%
%   *param* contains the following fields
%
%   * *param.order* : Degree of the Chebyshev approximation (default=100 
%     to yield a good approximation of the default Green's kernel). This
%     parameter is passed to gsp_filter_analysis.
%   * *param.regularize_epsilon*  : The regularized graph Laplacian is
%     $\bar{L}=L+\epsilon I$ (default = 0.005). A smaller epsilon may 
%     lead to better regularization, but will also require a higher order 
%     Chebyshev approximation. 
%
%   See also: gsp_filter_analysis 
% 
%   References: pesenson2009variational

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
    