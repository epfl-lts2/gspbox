function [sol, info] = gsp_prox_tv(x,gamma,G,param)
%GSP_PROX_TV Proximal TV operator for graphs signal
%   Usage:  sol = gsp_prox_tv(x, gamma, G, param)
%           sol = gsp_prox_tv(x, gamma, G)
%           [sol, info] = gsp_prox_tv(...)
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
%   This function computes the TV proximal operator for graphs. The TV norm
%   is the one norm of the gradient. The gradient is defined in the
%   function GSP_GRAD.
%
%   This function require the UNLocBoX to be executed.
%
%   GSP_PROX_TV(y, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||x||_TV
%
%   param is a Matlab structure containing the following fields:
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%
%     where  n(t) = f(x)+ 0.5 X-Z_2^2 is the objective function at iteration t*
%     by default, tol=10e-4.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.A : Forward operator (default: Id). This parameter allows to
%     solve the following problem
%
%        sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A x||_TV
%
%    param.At : Adjoint operator (default: Id).
%
%    param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%        ` ||A x||^2 <= nu * ||x||^2 
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.use_martrix : 1 use the matrix operation for the gradient.
%     This faster but requires more memory (default 1).
%
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used 
%
%   Demo: gsp_demo_graph_tv
%
%   See also: gsp_norm_tv gsp_grad
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/prox/gsp_prox_tv.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781



% Author: Nathanael Perraudin
% Date: 25 March 2014
% Testing: test_gsp_prox




if nargin < 3
    error('GSP_PROX_TV: You need to provide a graph!');
end

if nargin < 4, param=struct; end
if ~isfield(param, 'A'), param.A = @(x) x; end
if ~isfield(param, 'At'), param.At = @(x) x; end
if ~isfield(param, 'nu'), param.nu = 1; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'use_matrix'), param.use_matrix = 1; end;


% setting the function ftv
if ~isfield(G,'v_in')
    G = gsp_adj2vec(G);
    warning(['GSP_PROX_TV: To be more efficient you should run: ',...
        'G = gsp_adj2vec(G); before using this proximal operator.']);
end


if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_PROX_TV: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

param_l1.tight = 0;
param_l1.maxit = param.maxit;
param_l1.nu = 2*G.lmax*param.nu;

if param.use_matrix
    D = gsp_grad_mat(G);
    param_l1.A = @(x) D*param.A(x);
    param_l1.At = @(x) param.At(D'*x);
else
    param_l1.A = @(x) gsp_grad(G,param.A(x));
    param_l1.At = @(x) param.At(gsp_div(G,x));
end
param_l1.verbose = param.verbose;
param_l1.tol = param.tol;



[sol,info] = prox_l1(x,gamma,param_l1);

info.algo=mfilename;



end







