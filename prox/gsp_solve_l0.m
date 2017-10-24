function  [alpha, info]  = gsp_solve_l0(G, W, s, lambda, param )
%GSP_SOLVE_L0 Solve a BP problem using l0
%   Usage:  alpha = gsp_solve_l0(G, W, s, lambda );
%           alpha = gsp_solve_l0(G, W, s, lambda, param );
% 
%   Input parameters:
%       G       : Graph
%       W       : Filterbank
%       s       : Signal
%       lambda  : regularization parameter
%       param   : Optional structure of parameters
%   Ouptut parameters:
%       alpha   : Filterbank coefficients
%
%   This function solve the following minimization problem
%
%   .. argmin_alpha  lambda || alpha ||_0 + ||  W'alpha - s ||_2^2
%
%   .. argmin_\alpha \lambda \| \alpha \|_0  + \| W' \alpha - s \|_2^2
%
%   In *param*, you can optionaly set:
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   *info* is a Matlab structure containing optimization information from
%   the UNLocBoX
%

% Author: Nathanael Perraudin
% Date  : 8 February 2015
% Testing : test_gsp_solve_l0


if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end




Nf = length(W);
Ns = size(s,2);

% L0 minimization

paraml0.verbose = param.verbose -1;
fl0.prox = @(x,T) prox_l0(x,lambda * T,paraml0);
fl0.eval = @(x) lambda *sum(sum(abs(x)>0));


A = @(x) gsp_filter_synthesis(G,W,x,param);
At = @(x) gsp_filter_analysis(G,W,x,param);

[~, bound] = gsp_filterbank_bounds(G,W);

ffid.grad = @(x) 2* At(A(x)-s);
ffid.eval = @(x) norm(A(x)-s,'fro')^2;
ffid.beta = 2*bound;

% Starting point
if ~isfield(param,'guess'), param.guess = zeros(G.N*Nf,Ns); end

% Solver
paramsolver = param;

[alpha, info] = forward_backward(param.guess,fl0,ffid,paramsolver);




end

