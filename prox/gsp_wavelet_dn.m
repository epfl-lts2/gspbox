function [sol, info] = gsp_wavelet_dn(G, w, x, lambda, param)
%GSP_WAVELET_DN Wavelet denoising
%   Usage:  sol = gsp_wavelet_dn(G, w, x, lambda, param)
%           sol = gsp_wavelet_dn(G, w, x, lambda)
%           [sol, info] = gsp_wavelet_dn(...)
%
%   Input parameters:
%         G     : Graph structure
%         w     : Wavelet filterbank
%         x     : Signal to be denoised
%         lambda: Regularization parameter
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   This function will denoise a signal by solving the following convex
%   promblem:
%
%   .. sol = argmin_{z} || W^* z - x' ||_2^2 + lambda * ||z||_1
%
%   .. math::  sol = \min_{z} \| W^* z - x'\|_2^2 + \gamma  \|z\|_{TV}
%
%   Where $W$ is the frame associated to the filterbank $w$, $x'$ a part of
%   the signal to be denoised and *z* the wavelet coefficient. 
%
%   $x'$ consists of the high frequency part of *x*. It is obtained by
%   setting down to zero the low pass filter of the filerbank *w*.
%
%   This function require the UNLocBoX to be executed. You can download it
%   at http://unlcobox.sourceforge.net
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
%   * *param.tight* : 1 $W^*$ are both tight frame or 0 if not
%     (default = 0) 
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   * *param.method* : Solver to be used ('FISTA', 'ISTA', 'DG') By default
%     it is 'FISTA'. ('DG' is Douglas Rachford) 
%
%   * *param.gamma* : stepsize for the 'DG' algorithm
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
%   Demo: gsp_demo_wavelet_dn
%
%   See also: gsp_prox_tv gsp_solve_l1 gsp_solve_l0
%

% Author: Nathanael Perraudin
% Date: 25 March 2014
% Testing: 




if nargin<5
    param = struct;
end

if ~isfield(param, 'verbose'), param.verbose = 2; end
if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'method'), param.method = 'FISTA'; end
if ~isfield(param, 'gamma'), param.gamma = 1; end
if ~isfield(param, 'normalize'), param.normalize = 0; end

if param.normalize
    ntig = gsp_norm_tig(G,w(2:end));
    ntig = gsp_mat2vec(ntig);
else
    ntig = 1;
end

param_solver.maxit = param.maxit;
param_solver.verbose = param.verbose;
param_solver.tol = param.tol;

% % First method to handle low frequencies
% N = G.N;
% % Remove the low frequency part
% % 1) Analysis
% y_noise_c = gsp_filter_analysis(G,w,x);
% % % 2) Save low frequencies
% % y_noise_c_bf = zeros(size(y_noise_c));
% % y_noise_c_bf(1:N) = y_noise_c(1:N);
% % 3) Set low frequencies to 0
% y_noise_c(1:N) = 0;
% % 4) Inverse to find the signal without low frequencies
% y_noise_2 = gsp_filter_inverse(G,w,y_noise_c);
% y_noise_bf = x - y_noise_2;


% Remove the low frequency part
lf = gsp_design_expwin(G,0.3);
% 1) Analysis
y_noise_bf = gsp_filter_analysis(G,lf,x);
% 2) Save low frequencies
y_noise_2 = x-y_noise_bf;

% Sparsity of the coefficients -- prox l1
param_l1.verbose = param.verbose -1;
fc1.prox = @(x,T) prox_l1(x,T*lambda,param_l1);
fc1.eval = @(x) lambda*norm(x(:),1);
% Sparsity on the coefficients side -- prox l2
[~, B] = gsp_filterbank_bounds(G,w);

param_l2f.verbose = param.verbose - 1;
param_l2f.y = y_noise_2;
param_l2f.nu = 1/(2*B^2);
param_l2f.tight = param.tight;
param_l2f.ntig = ntig;
Wh = w(2:end); % drop the low pass filter
fl2.prox = @(x,T) gsp_prox_l2_filterbank(x, T, G, Wh, param_l2f);
fl2.eval = @(x) gsp_norm_l2_filterbank(G,Wh,x,param_l2f);
Wop = @(x) gsp_filter_analysis(G, Wh, x)./ntig;
Wtop = @(x) gsp_filter_synthesis(G, Wh, x./ntig);
fl2.grad = @(x) 2*Wop(Wtop(x)-y_noise_2);
fl2.beta = (2*B^2);

% starting point
x0 = Wop(x);

switch param.method
    case 'FISTA'
        param_solver.gamma = min(1/(2*B^2),0.5)*(min(ntig)^2);
        [denoise_y_wav2_c,info] = forward_backward(x0,fc1,fl2,param_solver);
    case 'ISTA'
        param_solver.gamma = min(1/(2*B^2),0.5)*(min(ntig)^2);
        param_solver.method = 'ISTA';
        [denoise_y_wav2_c,info] = forward_backward(x0,fc1,fl2,param_solver);
    case 'DG'
        param_solver.gamma = param.gamma;
        [denoise_y_wav2_c,info] = douglas_rachford(x0,fc1,fl2,param_solver);
    otherwise
        error('GSP_WAVELET_DN: Unknown method! ');
end


sol =Wtop(denoise_y_wav2_c)+y_noise_bf;
    
end
