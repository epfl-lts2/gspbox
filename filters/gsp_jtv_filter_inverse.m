function [s] = gsp_jtv_filter_inverse(G, filter, filtertype,c, param)
%GSP_JTV_FILTER_INVERSE Inverse operator of a joint timve_vertex filterbank
%   Usage:  s = gsp_jtv_filter_inverse(G, filter, c);
%           s = gsp_jtv_filter_inverse(G, filter, c, param);
%
%   Input parameters:
%         G          : Time-vertex Graph structure.
%         filter     : Cell array of time-vertex filters.
%         filtertype : Filter domain (ts,js,ts-array,js-array)
%         c          : Transform coefficients
%         param      : Optional parameter
%   Output parameters:
%         signal     : sythesis signal
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_filter_inverse.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

% Author: Francesco Grassi
% Testing: test_jtv_filter
% Date: September 2016

if nargin<5
    param = struct;
end

[dual_filter,filtertype] = gsp_jtv_design_can_dual(filter,filtertype);

s = gsp_jtv_filter_synthesis(G,dual_filter,filtertype,c,param);


end




