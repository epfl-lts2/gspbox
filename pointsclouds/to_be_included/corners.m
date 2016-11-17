function data = corners(N, scale, gapwidth, cornerwidth)

    if nargin < 1
        N = 1000;
    end
    if mod(N,8) ~= 0
        N = round(N/8) * 8;
    end

    if nargin < 2
        scale = 10;
    end
    if nargin < 3
        gapwidth = 2;
    end   
    if nargin < 4
        cornerwidth = 2;
    end

    perCorner = N/4;

    xplusmin = [ones(perCorner,1); -1*ones(perCorner,1); ones(perCorner,1); -1*ones(perCorner,1)];
    yplusmin = [ones(perCorner,1); -1*ones(2*perCorner,1); ones(perCorner,1)];
    
    horizontal = [xplusmin(1:2:end) * gapwidth + xplusmin(1:2:end) * scale .* rand(N/2,1), ...
                  yplusmin(1:2:end) * gapwidth + cornerwidth * yplusmin(1:2:end) .* rand(N/2,1), ...
                  floor((0:N/2-1)'/(perCorner*.5))];
       
    vertical = [xplusmin(2:2:end) * gapwidth + cornerwidth * xplusmin(2:2:end) .* rand(N/2,1), ...
                yplusmin(2:2:end) * gapwidth + yplusmin(2:2:end) * scale .* rand(N/2,1), ...
                floor((0:N/2-1)'/(perCorner*.5))];
    
    data=  [horizontal; vertical];

end
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/corners.php

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
end
