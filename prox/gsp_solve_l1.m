function  [alpha, info]  = gsp_solve_l1(G, W, s, lambda, param )
%GSP_SOLVE_L1 Solve a BP problem using l1
%   Usage:  alpha = gsp_solve_l1(G, W, s, lambda );
%           alpha = gsp_solve_l1(G, W, s, lambda, param );
% 
%   Input parameters:
%       G       : Graph
%       W       : Filterbank
%       s       : Signal
%       lambda  : Regularization parameter
%       param   : Optional structure of parameters
%   Ouptut parameters:
%       alpha   : Filterbank coefficients
%
%   This function solve the following minimization problem
%
%      argmin_alpha lambda || alpha ||_1  + || W'alpha - s ||_2^2
%
%      argmin_\alpha \lambda \| \alpha \|_1  + || W' \alpha - s ||_2^2
%
%   In param, you can optionaly set:
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%    param.normalize*: Use weight on the l1 norm to remove the effect of
%     the size of the atoms of W. Default 0
%    param.guess*: Initial guess for the starting point. (Default zeros)
%
%   info is a Matlab structure containing optimization information from
%   the UNLocBoX
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/prox/gsp_solve_l1.php

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
% Date  : 22 December 2014
% Testing : test_gsp_solve_l1


%% Todo, use primal dual solver!!!

if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'normalize'), param.normalize = 0; end








if param.normalize
    ntig = gsp_norm_tig(G,W);
    ntig = gsp_mat2vec(ntig);
    if iscell(ntig)
        ntig = cell2mat(ntig);
    end
else
    ntig = 1;
end
paraml1.weight = ntig;

paraml1.verbose = param.verbose -1;
paramsolver = param;

if lambda
    fl1.prox = @(x,T) prox_l1(x,lambda *T,paraml1);
    if param.normalize
        fl1.eval = @(x) lambda*sum(sum(abs(ntig.*x)));
    else
        fl1.eval = @(x) lambda*sum(sum(abs(x)));
    end
else
    fl1.prox = @(x,T) prox_l1(x,T,paraml1);
    if param.normalize
        fl1.eval = @(x) sum(sum(abs(ntig.*x)));
    else
        fl1.eval = @(x) sum(sum(abs(x)));
    end
    paramsolver = param;
    if ~isfield(paramsolver,'gamma'), paramsolver.gamma = 1/sqrt(G.N); end
    
end




if lambda
    A = @(x) gsp_filter_synthesis(G,W,x,param);
    At = @(x) gsp_filter_analysis(G,W,x,param);
    [~, bound] = gsp_filterbank_bounds(G,W);
    ffid.grad = @(x) 2* At(A(x)-s);
    ffid.eval = @(x) norm(A(x)-s,'fro')^2;
    ffid.beta = 2*bound;

else

    % Projection
    ffid.prox = @(x,T) gsp_proj_filterbank(x, T, G, W, s);
	ffid.eval = @(x) eps;

end

% Starting point
if ~isfield(param,'guess')
    N = G.N;
    Nf = numel(W);
    Ns = size(s,2);
    param.guess = zeros(N*Nf,Ns); 
end

paramsolver.stopping_criterion = 'rel_norm_primal';

% Solver
[alpha, info] = solvep(param.guess,{fl1,ffid},paramsolver);

end



