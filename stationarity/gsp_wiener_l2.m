function [sol, infos] = gsp_wiener_l2(G,y, A, At, psd, psd_noise, param)
%GSP_WIENER_l2 Solve wiener optimization problem with l2 fidelity term
%   Usage:  sol = gsp_wiener_l2(G, y, ffid, psd, psd_noise)
%           sol = gsp_wiener_l2(G, y, ffid, psd, psd_noise, param)
%           [sol, infos] = gsp_wiener_l2(...)
%
%   Input parameters:
%         G          : Graph (GSP structure)
%         x0         : Measurements (column vector)
%         A          : Operator (anonymous function)
%         At         : Adjoint operator (anonymous function)
%         psd        : PSD filter (anonymous function)
%         psd_noise  : PSD filter of the noise or single number
%         param      : Optional optimization parameters
%   Output parameters:
%         sol        : Solution
%         infos      : Convergence informations
%
%   This function solves the following wiener optimization problem:
%
%     .. argmin_x || A x - y ||_2^2 + || w(L) x ||_2^2 
%
%     .. math:: arg\min_x \| A x - y \|_2^2 + \| w(L) x \|_2^2 
%
%   Please refer to the reference for more information about this problem.
%   This function requires the UNLocBox to work.
%
%   Please refer to the function gsp_filter_analysis and solvep to know how
%   *param* can be set.
%
%   * *param.nu* : bound on the norm of the operator A (default: 1), i.e.
%
%     .. ` ||A x||^2 <= nu * ||x||^2 
%
%     .. math::  \|A x\|^2 \leq \nu  \|x\|^2   
%
%   References: perraudin2016stationary

% Author : Nathanael Perraudin
% Date: 6 January 2016


if nargin<7
    param = struct;
end

if ~isfield(param,'nu'), param.nu = 1; end
if ~isfield(param,'verbose'), param.verbose = 1; end

if isnumeric(psd_noise) && sum(abs(psd_noise))==0
    paramproj.A = A;
    paramproj.At = At;
    paramproj.epsilon = 1e-10;
    paramproj.y = y;
    paramproj.tight = 0;
    paramproj.verbose = param.verbose -1;
    ffid.prox = @(x,T) proj_b2(x,T,paramproj);
    ffid.eval = @(x) eps;

else
    % Fidelity term for Wiener optimization
    ffid.grad = @(x) 2*At(A(x)-y);
    ffid.eval = @(x) norm(A(x)-y,'fro')^2;
    ffid.beta = 2*param.nu;
end

[sol, infos] = gsp_wiener_optimization(G, y, ffid, psd, psd_noise, param);

end