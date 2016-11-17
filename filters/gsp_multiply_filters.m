function gm = gsp_multiply_filters(g1,g2)
%GSP_MULTIPLY_FILTERS Mutliply to filters
%   Usage: gm = gsp_multiply_filters(g1,g2);
%
%   Input parameters:
%       g1   : filterbank
%       g2   : filterbank
%
%   Output parameters:
%       gm  : multiplied filterbank
%
%   The resulting filter is gm(x) = g1(x) g2(x).
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_multiply_filters.php

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
% Date  : 30 September 2015

Nf1 = numel(g1);
Nf2 = numel(g2);

if ~iscell(g1)
    g1 = {g1};
end
if ~iscell(g2)
    g2 = {g2};
end

gm = cell(Nf1,Nf2);

for ii = 1:Nf1
    for jj = 1:Nf2
        gm{ii,jj} = @(x) g1{ii}(x).*g2{jj}(x);
    end
end


end
