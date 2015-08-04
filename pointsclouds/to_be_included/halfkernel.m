function data = halfkernel(N, minx, r1, r2, noise, ratio)

if nargin < 1
    N = 1000;
end
if mod(N,2) ~= 0
    N = N + 1;
end
if nargin < 2
    minx = -20;
end
if nargin < 3
    r1 = 20;
end
if nargin < 4
    r2 = 35;
end
if nargin < 5
    noise = 4;
end
if nargin < 6
    ratio = 0.6;
end

phi1 = rand(N/2,1) * pi;
inner = [minx + r1 * sin(phi1) - .5 * noise  + noise * rand(N/2,1) r1 * ratio * cos(phi1) - .5 * noise + noise * rand(N/2,1) ones(N/2,1)];
    
phi2 = rand(N/2,1) * pi;
outer = [minx + r2 * sin(phi2) - .5 * noise  + noise * rand(N/2,1) r2 * ratio * cos(phi2) - .5 * noise  + noise * rand(N/2,1) zeros(N/2,1)];

data = [inner; outer];
end
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/halfkernel.php

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
end
