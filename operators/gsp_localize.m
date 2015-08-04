function [ gt ] = gsp_localize(G, g, i)
%GSP_LOCALIZE Localize a kernel g to the node i
%   Usage: gt = gsp_localize(G, g, i);
%
%   Input parameters
%       G   : Graph
%       g   : kernel (or filterbank)
%       i   : Indices of vertex (int)
%   Output parameters
%       gt  : translate signal
%
%   This function localize the kernel g onto the node i. If g*
%   is a cell array, the localization will be done to each filter.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_localize.php

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
% Date  : 28 July 2014

f = zeros(G.N,1);
f(i) = 1;
gt = sqrt(G.N)*gsp_filter_analysis(G,g,f);

end


