function n = gsp_norm_l1_filterbank(G, W, x)
%GSP_NORM_L1_FILTERBANK Compute the l2 norm of the analysis coefficients
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_norm_l1_filterbank.php

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

%   Usage: n = gsp_norm_l1_filterbank(G, W, x);
%
%   Input parameters:
%       G   : Graph structure
%       W   : Filterbank (cell array of functions)
%       x   : coefficients
%       param: structure of optional paramter
%   Output parameters:
%       n   : L1 norm
%
%   `gsp_norm_l1_filterbank(G, W, x)` computes:
%
%   .. n = || W^* x ||_1
%
%   .. math::  n =  \|  W^* x \|_1
%
%

% Author: Nathanael Perraudin
% Date: 26 March 2014

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
warning(['GSP_NORM_L1_FILTERBANK: To be more efficient you should run: ',...
    'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

n = sum(sum(abs(gsp_filter_analysis(G,W,x))));

end

 

