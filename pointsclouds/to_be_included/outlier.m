function data = outlier(N, r, dist, outliers, noise)

    if nargin < 1
        N = 600;
    end
    if nargin < 2
        r = 20;
    end
    if nargin < 3
        dist = 30;
    end
    if nargin < 4
        outliers = 0.04;
    end
    if nargin < 5
        noise = 5;
    end

    N1 = round(N * (.5-outliers));
    N2 = N1;
    N3 = round(N * outliers);
    N4 = N-N1-N2-N3;

    phi1 = rand(N1,1) * pi;
    r1 = sqrt(rand(N1,1))*r;
    P1 = [-dist + r1.*sin(phi1) r1.*cos(phi1) zeros(N1,1)];

    phi2 = rand(N2,1) * pi;
    r2 = sqrt(rand(N2,1))*r;
    P2 = [dist - r2.*sin(phi2) r2.*cos(phi2) 3*ones(N2,1)];    
    
    P3 = [rand(N3,1)*noise dist+rand(N3,1)*noise 2*ones(N3,1)];    
    
    P4 = [rand(N4,1)*noise -dist+rand(N4,1)*noise ones(N4,1)];
    
    data = [P1; P2; P3; P4];

end
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/outlier.php

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
