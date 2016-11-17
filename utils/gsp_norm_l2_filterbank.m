function n = gsp_norm_l2_filterbank(G, W, x, param)
%GSP_NORM_L2_FILTERBANK Compute the l2 norm of the analysis coefficients
%   Usage: n = gsp_norm_l2_filterbank(G, W, x);
%
%   Input parameters:
%       G   : Graph structure
%       W   : Filterbank (cell array of functions)
%       x   : coefficients
%       param: structure of optional paramter
%   Output parameters:
%       n   : L2 norm
%
%   GSP_NORM_L2_FILTERBANK(G, W, x, param) computes:
%
%      n = || A W^* x -y ||_2^2
%
%   param is a Matlab structure containing the following fields:
%
%    param.A : Forward operator (default: Id).
%
%    param.y : measurements (default: 0).
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_norm_l2_filterbank.php

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
% Date: 26 March 2014

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
warning(['GSP_NORM_L1_FILTERBANK: To be more efficient you should run: ',...
    'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

if nargin < 4, param=struct; end

if ~isfield(param, 'y'), param.nu = 0; end
if ~isfield(param, 'A'), param.A = @(x) x; end

n = norm(param.A(gsp_filter_synthesis(G,W,x))-param.y,'fro')^2;

end

