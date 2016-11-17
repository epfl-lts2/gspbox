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
%        argmin_x || A x - y ||_2^2 + || w(L) x ||_2^2 
%
%   Please refer to the reference for more information about this problem.
%   This function requires the UNLocBox to work.
%
%   Please refer to the function gsp_filter_analysis and solvep to know how
%   param can be set.
%
%    param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%        ` ||A x||^2 <= nu * ||x||^2 
%
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/stationarity/gsp_wiener_l2.php

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
