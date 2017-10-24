function [sol, info] = gsp_proj_b2_filterbank(x, gamma, G, W, y, param)
%GSP_PROJ_B2_FILTERBANK Projection on the B2 ball for a filterbank
%   Usage:  sol = gsp_proj_b2_filterbank(x, T, G, W, param);
%           sol = gsp_proj_b2_filterbank(x, T, G, W);
%           [sol, info] = gsp_proj_b2_filterbank(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Compatibility parameter
%         G     : Graph structure
%         W     : Filterbank (cell array of functions)
%         y     : Measurements
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   This function require the UNLocBoX to be executed.
%
%   `gsp_proj_b2_filterbank(x, gamma, G, W, param)` can solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2  such that || W A x -y ||_2 < epsilon 
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 \text{ s. t. }  \| W A x - y\|_2 < \epsilon 
%
%   Where $W$ is the linear analysis operator associated with the
%   filterbank. In this case, we solve the problem on the signal side.
%   Alternatively, we can also solve it on the coefficient side.
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2  such that || A W^* x -y ||_2 < epsilon 
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 \text{ s. t. }  \| A W^* x - y\|_2 < \epsilon 
%
%   You can select the problem you want to solve by setting *param.type* to
%   'coefficient' or 'signal'.
%
%   param is a Matlab structure containing the following fields:
%   
%   * *param.tight* : 1 if $A W^*$ is a tight frame or 0 if not (default = 0)
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
%   * *param.epsilon* : Radius of the L2 ball (default = 1e-3).
%
%   * *param.type* : 'coefficient' or signal 'signal' select the problem
%     type. (default 'signal').
%
%   * *param.normalize*: normalize the frame (for testing purpose only)
%     this is not efficient
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
%   See also:  gsp_prox_l1_filterbank gsp_prox_l2_filterbank gsp_prox_tv
%


% Author: Nathanael Perraudin
% Date: 25 March 2014
%

if param.epsilon == 0
    warning('This function could be greatly accelerated using dual filters, please ask Nathanael.')
end



if nargin < 5
    error('GSP_PROJ_B2_FILTERBANK: Not enought input arguments');
end

if nargin < 6, param=struct; end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'At'), param.At = param.A; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'type'), param.type = 'signal'; end
if ~isfield(param, 'normalize'), param.normalize = 0; end

if iscell(G)
    Gtemp = G{1};
else
    Gtemp = G;
end

if ~isfield(Gtemp,'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_PROJ_B2_FILTERBANK: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

if param.normalize
    ntig = gsp_norm_tig(G,W);
    ntig = gsp_mat2vec(ntig);
    if iscell(ntig)
        ntig = cell2mat(ntig);
    end
else
    ntig = 1;
end

[~, B] = gsp_filterbank_bounds(G,W);
if iscell(B)
    B = sum(cell2mat(B));
end
param_b2.epsilon = param.epsilon;
param_b2.tight = param.tight;
param_b2.maxit = param.maxit;
param_b2.nu =  B*param.nu/min(ntig)*sqrt(2);
param_b2.y = y;
switch param.type
    case 'coefficient'
        param_b2.A= @(x) syn_c(G,W,x,ntig,param);
        param_b2.At = @(x) ana_c(G,W,x,ntig,param);
    case 'signal'
        param_b2.A= @(x) ana_s(G,W,x,ntig,param);
        param_b2.At = @(x) syn_s(G,W,x,ntig,param);
    otherwise
        error('GSP_PROJ_B2_FILTERBANK: Unknown type of problem!')
end
param_b2.verbose = param.verbose;


[sol,info] = proj_b2(x,gamma,param_b2);

info.algo=mfilename;


end


% Here I handle operator with cell array. It is a pain in the ass
function s = syn_c(G,W,x,ntig,param)
    fon = @(G, W, x ,ntig) param.A(gsp_filter_synthesis(G,W,x./ntig));
    s = func(fon, G, W, x, ntig, 1);
end

function s = ana_c(G,W,x,ntig,param)
    fon = @(G, W, x ,ntig) gsp_filter_analysis(G,W,param.At(x))./ntig;
    s = func(fon, G, W, x, ntig, 0);
end

function s = syn_s(G,W,x,ntig,param)
    fon = @(G, W, x ,ntig) param.At(gsp_filter_synthesis(G,W,x./ntig));
    s = func(fon, G, W, x, ntig, 1);
end

function s = ana_s(G,W,x,ntig,param)
    fon = @(G, W, x ,ntig) gsp_filter_analysis(G,W,param.A(x))./ntig;
    s = func(fon, G, W, x, ntig, 0);
end

function s = func(fon, G, W, x, ntig, syn)
if iscell(G) % Cell case
    N = G{1}.N;
    NF = N * length(W{1});
    Nx = size(x,2);
    NG = numel(G);
    if syn
        s = zeros(N,Nx);
        for ii = 1:NG
            vec = (1:NF) + (ii-1) * NF;
            s = s + fon(G{ii},W{ii},x(vec,:),getntig(ntig, vec));
        end
    else
        s = zeros(NF*NG,Nx);
        for ii = 1:NG
            vec = (1:N) + (ii-1) * N;
            vec2 = (1:NF) + (ii-1) * NF;
            s(vec2,:) = fon(G{ii},W{ii},x,getntig(ntig, vec));
        end
    end
else % simple case
    s = fon(G,W,x,ntig);
end
end

function ntig = getntig(ntig,vec)
    if numel(ntig)>1
        ntig = ntig(vec);
    end
end

