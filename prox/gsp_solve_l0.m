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
%      argmin_alpha  lambda || alpha ||_0 + ||  W'alpha - s ||_2^2
%
%      argmin_\alpha \lambda \| \alpha \|_0  + \| W' \alpha - s \|_2^2
%
%   In param, you can optionaly set:
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%   info is a Matlab structure containing optimization information from
%   the UNLocBoX
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/prox/gsp_solve_l0.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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


