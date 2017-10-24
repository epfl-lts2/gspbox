function [sol, infos] = gsp_ml_rls(xl, y, k, tau, A, At, param)
%GSP_ML_RLS Manifold Learning regularized least square
%   Usage: sol = gsp_ml_rls(xl, y, k, tau);
%          sol = gsp_ml_rls(xl, y, k, tau, A, At);
%          sol = gsp_ml_rls(xl, y, k, tau,  A, At, param);
%          [sol, infos] = gsp_ml_rls(...)
%   
%   Input parameters:
%       xl      : labeled points
%       y       : labels
%       k       : kernel
%       tau     : regularization parameters
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
%   ..  argmin_alpha  || A (K alpha) - y ||_2^2 + tau *alpha^T K alpha
%
%   If tau is set to zero, then the following problem is solved
%
%   ..  argmin_alpha alpha^T K alpha s. t.  A (K alpha) = y
%

% Author: Nathanael Perraudin
% Date  : 8 decembre 2014


if nargin<5
    A = @(x) x;
end


if nargin<6
    At = A;
end

if nargin<7
    param = struct;
end


if ~isfield(param, 'tol'), param.tol = 1e-6; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
    
    

%N = size(xl,2);

% Evaluate the kernel on the data points
K = gsp_rkhs_evaluate(k,xl);
nu = norm(K);

[N,M] = size(y);
Nk = length(k);
alpha_in = zeros(N*Nk,M);

if tau >0
    
    fp.eval = @(x) tau * sum(norm_rkhs( K,x ));
    fp.grad = @(x) tau * grad_rkhs( K,x);
    fp.beta =   2*nu*tau;

    ffid.eval = @(x) norm(A(K*x)-y,'fro')^2;
    ffid.grad = @(x) 2*K'*At(A(K*x)-y);
    ffid.beta = 2*nu^2*param.nu^2;


    
    
%     paramfid.A = @(x) A(K*x);
%     paramfid.At = @(x) K'*At(x);
%     paramfid.nu = nu^2;
%     paramfid.tight = 0;
%     paramfid.y = y;
%     paramfid.verbose = param.verbose -1;
%     ffid.eval = @(x) norm(A(K*x)-y,'fro')^2;
%     ffid.prox = @(x,T) prox_l2(x,T,paramfid);
    
    [sol,infos] = solvep(alpha_in, {fp,ffid}, param);

%     param.gamma = 0.5/(tau*nu);
%     [sol,infos] = forward_backward(alpha_in, ffid, fp, param);
else
    
    fp.eval = @(x) sum(norm_rkhs( K,x ));
    fp.grad = @(x) grad_rkhs( K,x);
    fp.beta =   2*nu;

    
    paramproj.A = @(x) A(K*x);
    paramproj.At = @(x) K'*At(x);
    paramproj.nu = nu^2*param.nu^2;
    paramproj.tight = 0;
    paramproj.verbose = param.verbose-1;
    paramproj.y = y;
    paramproj.maxit = 50;
    ffid.eval = @(x) eps;
    ffid.prox = @(x,T) proj_b2(x,T,paramproj);


    [sol,infos] = forward_backward(alpha_in, ffid, fp, param);

end




end