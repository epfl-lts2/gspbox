function [ gt ] = gsp_localize(G, g, n,param)
%GSP_LOCALIZE Localize a kernel g to the node n
%   Usage: gt = gsp_localize(G, g, n);
%
%   Input parameters
%       G   : Graph
%       g   : kernel (or filterbank)
%       n   : Indices of vertex (int)
%       param: Optional parameters
%   Output parameters
%       gt  : translate signal
%
%   This function localize the kernel g onto the node i. If g*
%   is a cell array, the localization will be done to each filter.
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/operators/gsp_localize.html

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

% Author: Nathanael Perraudin
% Date  : 28 July 2014

if nargin <4
    param = struct;
end

f = zeros(G.N,numel(n));
for ii = 1:numel(n)
    f(n(ii),ii) = 1;
end
gt = gsp_filter_analysis(G,g,f,param);

end


