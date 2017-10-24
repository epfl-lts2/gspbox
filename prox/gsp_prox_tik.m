function [sol, info] = gsp_prox_tik(x,gamma,G,param)
%GSP_PROX_TIK Proximal tikhonov operator for graphs
%   Usage:  sol = gsp_prox_tik(x, gamma, G, param)
%           sol = gsp_prox_tik(x, gamma, G)
%           [sol, info] = gsp_prox_tik(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         G     : Graph structure
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   This function require the UNLocBoX to be executed.
%
%   `gsp_prox_tik(y, gamma, param)` solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * || nabla x||_2^2
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \| \nabla x\|_2^2
%
%   Note the nice following relationship
%
%   ..   x' L x = || nabla x||_2^2
%
%   .. math::  x^T L x = || \nabla x||_2^2
%
%   param is a Matlab structure containing the following fields:
%   
%   * *param.tol* : is stop criterion for the loop. The algorithm stops if
%
%     ..  (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     .. math:: \frac{  n(t) - n(t-1) }{ n(t)} < tol,
%
%     where  $n(t) = f(x)+ 0.5 \|x-z\|_2^2$ is the objective function at iteration *t*
%     by default, `tol=10e-4`.
%
%   * *param.maxit* : max. nb. of iterations (default: 200).
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.A* : Forward operator (default: Id). This parameter allows to
%     solve the following problem
%
%     .. sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * || nabla A x||_2^2
%
%     .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 + \gamma  \|\nabla A x\|_2^2 
%
%   * *param.At* : Adjoint operator (default: Id).
%
%   * *param.pcg* : Use the fast PCG algorithm (default 1).
%
%   * *param.use_martrix* : 1 use the matrix operation for the gradient.
%     This faster but requires more memory (default 1).
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
%   See also:  gsp_prox_tv
%



% Author: Nathanael Perraudin
% Date: 25 March 2014
% Testing: test_gsp_prox


% Start the time counter
t1 = tic;

if nargin < 3
    error('GSP_PROX_TIK: You need to provide a graph!');
end

if nargin < 4, param=struct; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'pcg'), param.pcg = 1; end
if ~isfield(param, 'use_matrix'), param.use_matrix = 1; end;
if ~isfield(param, 'order'), param.order = 30; end;



if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    if param.verbose
    warning(['GSP_PROX_TV: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
    end
end

if ~isfield(param,'A')


    h =@(x) 1./(1+2*gamma*x);
    paramfilter.order = param.order;
    sol = gsp_filter_analysis(G,h,x,paramfilter);

    obj = sum(gsp_norm_tik(G,sol)) + 0.5*norm(x-sol)^2;

    if param.verbose
        fprintf(['  GSP_Prox_tik: 0.5*||x - z||_2^2 +', ...
                ' gamma * || nabla x||_2^2 = %e\n'], obj);
    end
    
    if isa(x,'single')
        sol = single(sol); 
    end
    
    info.algo=mfilename;
    info.iter=1;
    info.final_eval = obj;
    info.crit='--';
    info.time=toc(t1);
else

    if ~isfield(param,'At')
        error('GSP_PROX_TIK: Please provide param.At')
    end
    

    if ~isfield(G,'v_in')
        G = gsp_adj2vec(G);
        warning(['GSP_PROX_TV: To be more efficient you should run: ',...
            'G = gsp_adj2vec(G); before using this proximal operator.']);
    end

    param_l2.tight = 0;
    param_l2.maxit = param.maxit;
    param_l2.tol = param.tol;
    param_l2.nu = 2*G.lmax*param.nu;

    if param.use_matrix
        D = gsp_grad_mat(G);
        param_l2.A = @(x) D*param.A(x);
        param_l2.At = @(x) param.At(D'*x);
    else
        param_l2.A = @(x) gsp_grad(G,param.A(x));
        param_l2.At = @(x) param.At(gsp_div(G,x));
    end
    param_l2.verbose = param.verbose;
    param_l2.pcg = param.pcg;


    [sol,info] = prox_l2(x,gamma,param_l2);
    
    info.algo=mfilename;
 

end



end






