function [f_interpolated,varargout]=gsp_interpolate_old(f_subsampled,G,keep_inds,varargin)
%GSP_INTERPOLATE Interpolation of a graph signal
%   Usage:  f_interpolated=gsp_interpolate(f_subsampled,G,keep_inds);
%           f_interpolated=gsp_interpolate(f_subsampled,G,keep_inds,param);
%
%   Input parameters:
%         f_subsampled     : A graph signal on the graph G (with G.N
%                            vertices), with the values on $V_1^c$ set to 0
%         G                : Graph structure.
%         keep_inds        : The vertex indices of $V_1$; i.e., those 
%                            where the signal values are available to use 
%                            in the interpolation.
%   Output parameters:
%         f_interpolated   : Interpolated graph signal on G.
%   Additional parameters:
%         param.use_exact           : To use exact graph spectral filtering instead of the Chebyshev approximation.
%         param.order          : Degree of the Chebyshev approximation (default=30).
%         param.regularize_epsilon  : The regularized graph Laplacian is $\bar{L}=L+\epsilon I$. 
%                                     A smaller epsilon may lead to better regularization, but will also require a higher order Chebyshev approximation.
%
%   'gsp_interpolate(f_subsampled,G,keep_inds)' performs the interpolation:
%
%   .. f_interpolated = \Phi_{V_1} \Phi_{V_1,V_1}^{-1} f_{V_1} = \Phi_{V_1} [\bar{L}_{V_1,V_1}-\bar{L}_{V_1,V_1^c}(\bar{L}_{V_1^c,V_1^c})^{-1}\bar{L}_{V_1^c,V_1}] f_{V_1} 
%
%   .. math:: f_{interpolated}=\Phi_{{\cal V}_1} \Phi_{{\cal V}_1,{\cal V}_1}^{-1} f_{{\cal V}_1} = \Phi_{{\cal V}_1} \left[\bar{{\cal L}}_{{\cal V}_1,{\cal V}_1}-\bar{{\cal L}}_{{\cal V}_1,{\cal V}_1^c}(\bar{{\cal L}}_{{\cal V}_1^c,{\cal V}_1^c})^{-1}\bar{{\cal L}}_{{\cal V}_1^c,{\cal V}_1}\right] f_{{\cal V}_1}
%

%
%   See also:  
%
%   Demos:  
% 
%   References: I. Pesenson, "Variational splines and Paley-Wiener spaces
%   on combinatorial graphs," Constr. Approx., vol. 29, no. 1, pp.
%   31-21, Feb. 2009.

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
% Read input parameters
if nargin>3
    param=varargin{1};
else
    param=0;
end

if ~isfield(param,'use_exact');
    use_exact=( isfield(G,'U') && isfield(G,'e') );
else
    use_exact=param.use_exact;
end

elim_inds=setdiff(1:G.N,keep_inds);


if ~isfield(param, 'regularize_epsilon')
    epsilon=.005;
else
    epsilon = param.regularize_epsilon;
end

% Compute coefficients alpha for each of the Green's functions translated
% to center vertices in V_1
regularized_L= G.L+epsilon*eye(G.N);
green_kernel=@(x) 1./(x+epsilon);
alpha=regularized_L(keep_inds,keep_inds)*f_subsampled(keep_inds)-regularized_L(keep_inds,elim_inds)*(regularized_L(elim_inds,elim_inds)\(regularized_L(elim_inds,keep_inds)*f_subsampled(keep_inds)));

% Now compute f_interp=\Phi_{V_1} \alpha
if use_exact 
    if ( ~isfield(G,'e') || ~isfield(G,'U') )
        G=gsp_full_eigen(G);
        varargout{1}=G;
    end
    green_functions=G.U*diag(green_kernel(G.e))*G.U';
    f_interpolated=green_functions(:,keep_inds)*alpha;

else % use Chebyshev approximation method

    if ~isfield(param,'order');
        order=100;
    else
        order=param.order;
    end

    if ~isfield(G,'lmax');
        G.lmax=sgwt_rough_lmax(G.L);
    end

    cheb_coeffs=sgwt_cheby_coeff(green_kernel,order,order+1,[0,G.lmax]);
    alpha_upsampled=zeros(G.N,1);
    alpha_upsampled(keep_inds)=alpha;
    f_interpolated=sgwt_cheby_op(alpha_upsampled,G.L,cheb_coeffs,[0,G.lmax]);
end
    
