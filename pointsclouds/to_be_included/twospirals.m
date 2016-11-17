function data = twospirals(N, degrees, start, noise)
% Generate "two spirals" dataset with N instances.
% degrees controls the length of the spirals
% start determines how far from the origin the spirals start, in degrees
% noise displaces the instances from the spiral. 
%  0 is no noise, at 1 the spirals will start overlapping
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/twospirals.php

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

    if nargin < 1
        N = 2000;
    end
    if nargin < 2
        degrees = 570;
    end
    if nargin < 3
        start = 90;
    end
    if nargin < 5
        noise = 0.2;
    end  
    
    deg2rad = (2*pi)/360;
    start = start * deg2rad;

    N1 = floor(N/2);
    N2 = N-N1;
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;   
    d1 = [-cos(n).*n + rand(N1,1)*noise sin(n).*n+rand(N1,1)*noise zeros(N1,1)];
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;      
    d2 = [cos(n).*n+rand(N2,1)*noise -sin(n).*n+rand(N2,1)*noise ones(N2,1)];
    
    data = [d1;d2];
end
