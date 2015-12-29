function [x,y] = gsp_twospirals(N, degrees, start, noise)
%GSP_TWOSPIRALS Generate "two spirals" dataset with N instances
%   Usage:  gsp_twospirals(N, degrees, start, noise);
%           gsp_twospirals(N, degrees, start);
%           gsp_twospirals(N, degrees);
%           gsp_twospirals(N);
%           gsp_twospirals();
%
%   Input arguments:
%           N       : Number of points (default 2000)
%           degrees : The length of the spirals in degree (default 570)
%           start   : how far from the origin the spirals start, in degrees (default 90)
%           noise   : Noise level (default 0.2) 
%
%   Output arguments:
%           x       : Position of the points
%           y       : label
%           
%   Note that for noise=0, there is no noise and at noise=1 the spirals
%   will start overlapping 
%
%   This function is adaptated from:
%   http://stackoverflow.com/questions/16146599/create-artificial-data-in-matlab
%   and 
%   http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
%
%   Example:
%
%       [x,y] = gsp_twospirals();
%       figure(); 
%       plot(x(y>0,1),x(y>0,2),'xr');
%       hold on
%       plot(x(y<=0,1),x(y<=0,2),'ob');
%       hold off
%       
%   See also: 
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/gsp_twospirals.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Date: 15 August 2015
% Adaptated by: Nathanael Perraudin

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
    x = data(:,1:2);
    y = data(:,3);
end
