function [sol, infos] = gsp_ml_rls_tv(G, xl, y, k, tau,lambda, A, At, param)
%GSP_ML_RLS_TV Manifold Learning regularized least square with TV regularization
%   Usage: sol = gsp_ml_rls_tv(G, xl, y, k, tau,lambda);
%          sol = gsp_ml_rls_tv(G, xl, y, k, tau, lambda, A, At);
%          sol = gsp_ml_rls_tv(G, xl, y, k, tau, lambda, A, At, param);
%          [sol, infos] = gsp_ml_rls_tv(...)
%   
%   Input parameters:
%       G       : Graph
%       xl      : labeled points
%       y       : labels
%       k       : kernel
%       tau     : regularization parameters
%       lambda  : regularization parameters
%       A       : Operator
%       At      : Adoint operator
%       param   : Optional parameters
%   Output parameters:
%       sol     : solution of the problem (kernel coefficients)
%       infos   : convergence info
%       
%   *param* is a structure of optional argument given to the solver
%   gradient_descent. Please see the function gradient descent for more
%   information. 
%
%   In *param*, you also have to set an upperbound for the operator A as
%   param.nu!
%
%   This function solves the following problem:
%
%   ..  argmin_alpha  || A (K alpha) - y ||_2^2 
%                       + tau *alpha^T K alpha 
%                       + lambda || L K alpha ||_TVG
%
%   If tau is set to zero, then the following problem is solved
%
%   ..  argmin_alpha alpha^T K alpha
%                      + lambda || L K alpha ||_TVG
%                      s. t.  A (K alpha) = y
%

% Author: Nathanael Perraudin
% Date  : 8 decembre 2014


if nargin<7
    A = @(x) x;
end


if nargin<8
    At = A;
end

if nargin<9
    param = struct;
end


if ~isfield(param, 'tol'), param.tol = 1e-6; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
    
    

N = size(xl,2);

% Evaluate the kernel on the data points
K = gsp_rkhs_evaluate(k,xl);
nu = norm(K);

alpha_in = zeros(N,size(y,2));

if ~isfield(G,'D');
    G = gsp_adj2vec(G);
end
if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);
end
    
paramtv.verbose = param.verbose - 1;
ftv.eval = @(x) sum(lambda*gsp_norm_tv(G,K*x));
ftv.prox = @(x,T) gsp_prox_tv(x,lambda*T,G,paramtv);


if tau >0
    
    fp.eval = @(x) tau * sum(gsp_norm_tik(K,x));
    fp.grad = @(x) 2*tau*K*x;

    ffid.eval = @(x) norm(A(K*x)-y,'fro')^2;
    ffid.grad = @(x) 2*K'*At(A(K*x)-y);
    
    ftot.grad = @(x) ffid.grad(x) + fp.grad(x);
    ftot.eval = @(x) ffid.eval(x) + fp.eval(x);


    param.gamma = 0.5/(tau*nu+nu^2*param.nu^2+lambda*nu*G.lmax);

    [sol,infos] = forward_backward(alpha_in, ftv, ftot, param);
else
    
    fp.eval = @(x)  x'*K*x;
    fp.grad = @(x) 2*K*x;
    
    
    paramproj.A = @(x) A(K*x);
    paramproj.At = @(x) K'*At(x);
    paramproj.nu = nu^2*param.nu^2;
    paramproj.tight = 0;
    paramproj.verbose = param.verbose-1;
    paramproj.y = y;
    paramproj.maxit = 50;
    ffid.eval = @(x) eps;
    ffid.prox = @(x,T) proj_b2(x,T,paramproj);

    param.gamma = 0.5/(nu+lambda*nu*G.lmax);

    [sol,infos] = generalized_forward_backward(alpha_in, {ffid,ftv}, fp, param);

end




end