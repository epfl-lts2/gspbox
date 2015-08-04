function data = crescentfullmoon(N, r1, r2, r3)

if nargin < 1
    N = 1000;
end
if mod(N,4) ~= 0
    N = round(N/4) * 4;
end
if nargin < 2
    r1 = 5;
end
if nargin < 3
    r2 = 10;
end
if nargin < 4
    r3 = 15;
end

N1 = N/4;
N2 = N-N1;

phi1 = rand(N1,1) * 2 * pi;
R1 = sqrt(rand(N1, 1));
moon = [cos(phi1) .* R1 * r1 sin(phi1) .* R1 * r1 zeros(N1,1)];

d = r3 - r2;
phi2 = pi + rand(N2,1) * pi;
R2 = sqrt(rand(N2,1));
crescent = [cos(phi2) .* (r2 + R2 * d) sin(phi2) .* (r2 + R2 * d) ones(N2,1)];

data = [moon; crescent];
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/crescentfullmoon.php

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
data = [moon; crescent];
