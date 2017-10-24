function [sol, info] = gsp_prox_l2_filterbank(x, gamma, G, W, param)
%GSP_PROX_L2_FILTERBANK Proximal L2 operator for a filterbank
%   Usage:  sol = gsp_prox_l2_filterbank(x, T, G, W, param);
%           sol = gsp_prox_l2_filterbank(x, T, G, W);
%           [sol, info] = gsp_prox_l2_filterbank(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter
%         G     : Graph structure
%         W     : Filterbank (cell array of functions)
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function require the UNLocBoX to be executed.
%
%   `gsp_prox_l2_filterbank(x, gamma, G, W, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2  || A W^* x - y ||_2^2
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \| A W^* x - y\|_2^2 
%
%   Where $W$ is the linear analysis operator associated with the
%   filterbank.
%
%   param is a Matlab structure containing the following fields:
%   
%   * *param.tight* : 1 if A is a tight frame or 0 if not (default = 0)
%
%   * *param.y* : measurements (default: 0).
%   
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f(x)+ 0.5 \|x-z\|_2^2$ is the objective function at
%     iteration *t* by default, `tol=10e-4`.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.weights* : weights for a weighted L2-norm (default = 1)
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.At* : Adjoint operator (default: A).
%
%   * *param.nu* : bound on the norm of the operator A (default: 1), i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2 
%
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%
%   * *info.iter* : Number of iteration
%
%   * *info.time* : Time of exectution of the function in sec.
%
%   * *info.final_eval* : Final evaluation of the function
%
%   * *info.crit* : Stopping critterion used 
%
%
%   See also:  gsp_prox_l1_filterbank gsp_norm_l2_filterbank gsp_prox_tv
%



% Author: Nathanael Perraudin
% Date: 25 March 2014
%


if nargin < 3
    error('GSP_PROX_L2_FILTERBANK: You need to provide a graph!');
end

if nargin < 4, param=struct; end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'y'), param.nu = 0; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'At'), param.At = param.A; end
if ~isfield(param, 'ntig'), param.ntig = 1; end

% setting the function ftv

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_PROX_L2_FILTERBANK: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

[~, B] = gsp_filterbank_bounds(G,W);

param_l2.tight = param.tight;
param_l2.maxit = param.maxit;
param_l2.nu =  B*param.nu/min(param.ntig);
param_l2.y = param.y;
param_l2.A= @(x) param.A(gsp_filter_synthesis(G,W,x./param.ntig));
param_l2.At = @(x) gsp_filter_analysis(G,W,param.At(x))./param.ntig;
param_l2.verbose = param.verbose;


[sol,info] = prox_l2(x,gamma,param_l2);

info.algo=mfilename;



end






