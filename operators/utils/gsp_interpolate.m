function s_pred = gsp_interpolate(Gh,Gl,coeff,param)
%GSP_INTERPOLATE Interpolate lower coefficient
%   Usage: s_pred = gsp_interpolate(G,K_reg,ind,g,coeff);
%
%   Input parameters:
%       Gh      : Upper graph
%       Gl      : Lower graph
%       coeff   : Coefficients
%       param   : Optional parameters for gsp_filter_anlysis
%   Ouptut parameters:
%       s_pred  : Predicted signal
%
%   This function is for internal use only
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/utils/gsp_interpolate.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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
% Date   : 14 August 2014
% Testing: test_pyramid

if nargin < 4
    param = struct;
end
    
if ~isfield(param,'cheb_order'), param.cheb_order = 100; end

alpha = Gl.pyramid.K_reg * coeff;
s_pred = zeros(Gh.N,1);
s_pred(Gl.pyramid.ind) = alpha;
s_pred = gsp_filter_analysis(Gh,Gl.pyramid.green_kernel ,s_pred,param);
    
end
