function [sol, info] = gsp_proj_filterbank(x, ~ , G, W, y, param)
%GSP_PROJ_FILTERBANK Projection onto the synthesis coefficients
%   Usage:  sol = gsp_proj_filterbank(x, 0, G, W, y, param);
%           sol = gsp_proj_filterbank(x, 0, G, W, y);
%           [sol, info] = gsp_proj_filterbank(...)
%
%   Input parameters:
%         x     : Input signal
%         G     : Graph structure
%         W     : Filterbank (cell array of functions)
%         y     : Measurements
%         param : Structure of optional parameters
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   `gsp_proj_filterbank(x, gamma, G, W, param)` can solves:
%
%   .. sol = argmin_{z} 0.5*||x - z||_2^2  such that W^* x =y 
%
%   .. math::  sol = \min_{z} \frac{1}{2} \|x - z\|_2^2 \text{ s. t. }  W^* x = y 
%
%   Where $W$ is the linear analysis operator associated with the
%   filterbank. 
%
%   The function can use different techniques
%   
%   * 'exact' : if the Fourier basis is computed, go for this one
%   * 'cheby' : use the pseudo-inverse filters of the filterbank with
%     chebyshev approximation. It works well for well-conditionned
%     filterbanks.
%   * 'lanczos' : use the pseudo-inverse filters of the filterbank with
%     lanczos approximation. It works well for well-conditionned
%     filterbanks.
%   * 'proj_b2': scallable and robust way to do it. However, the
%     convergence maybe slow and might require a lot of filtering
%     operations.
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%   * *param.eps* : tolerance for the pseudo inverse method
%   * *param.proj_method*: selected method
%
%   info is a Matlab structure containing the following fields:
%
%   * *info.algo* : Algorithm used
%   * *info.iter* : Number of iteration
%   * *info.time* : Time of exectution of the function in sec.
%   * *info.final_eval* : Final evaluation of the function
%   * *info.crit* : Stopping critterion used 
%
%
%   See also:  gsp_solve_l1 gsp_proj_b2_filterbank 
%


% Author: Nathanael Perraudin
% Date: 25 March 2014
% Testing: test_gsp_proj_filerbank




if nargin < 5
    error('GSP_PROJ_FILTERBANK: Not enought input arguments');
end

if nargin < 6, param=struct; end

if ~isfield(param, 'eps'), param.eps = 1e-8; end
if ~isfield(param, 'proj_method'), 
    if gsp_check_fourier(G)
        param.proj_method = 'exact';
    else
        [A,B] = gsp_filterbank_bounds(G,W);
        if B/A > 20
            warning('Filterbank ill-conditioned, going for a slow method');
            param.proj_method = 'primal_dual';
        else
            param.proj_method = 'cheby';
        end
    end
end




t = tic;

if size(y,2)==1 && size(x,2)>1
    y = repmat(y,1,size(x,2));
end
if isfield(param,'method')
    param2 = rmfield(param,'method');
else 
    param2 = param;
end
if strcmp(param.proj_method, 'proj_b2') || strcmp(param.proj_method, 'primal_dual')
    [~,B] = gsp_filterbank_bounds(G,W);
    if ~isfield(param,'paramproj'), param.paramproj = struct; end
    paramproj = param.paramproj;
    paramproj.A = @(x) gsp_filter_synthesis(G,W,x,param2);
    paramproj.At = @(x) gsp_filter_analysis(G,W,x,param2);
    paramproj.nu = B^2;
    paramproj.epsilon = sqrt(G.N)*param.eps;
    paramproj.y = y;
    paramproj.method  = param.proj_method;

    [sol,info] = proj_linear_eq(x,0,paramproj);
    
else
    if strcmp(param.proj_method, 'cheby')  || strcmp(param.proj_method, 'lanczos')
        [A,B] = gsp_filterbank_bounds(G,W);
        if B/A > 20
            warning('Filterbank ill-conditioned, check your solution!');
        end
    end
    Wd = gsp_design_can_dual(W,param.eps);
    sol =  x - gsp_filter_analysis(G,Wd,(gsp_filter_synthesis(G,W,x,param2)-y),param2);
    
    info.iter = 1;
    info.final_eval = 0;
    info.crit = 'Direct computation';
end

    info.time = toc(t);
    info.algo = param.proj_method;


end


