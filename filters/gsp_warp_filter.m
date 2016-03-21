function gw = gsp_warp_filter(g,w)
%GSP_WARP_FILTER Warp the filterbank g with the filter w
%   Usage: gw = gsp_warp_filter(g,w);
%
%   Input parameters:
%       g   : filterbank
%       w   : warping filter
%
%   Output parameters:
%       gw  : warped filterbank
%
%   The resulting filter gw is gw(x)=w(g(x)).
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_warp_filter.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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
% Date  : 30 September 2015

Nf = numel(g);

if ~iscell(g)
    g = {g};
end

gw = cell(Nf,1);

for ii = 1:Nf
    gw{ii} = @(x) w(g{ii}(x));
end


end
