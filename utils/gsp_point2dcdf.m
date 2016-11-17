function [x,y] = gsp_point2dcdf(v)
%GSP_POINT2DCDF points to discrete continuous density function
%   Usage : [x,y] = gsp_point2dcdf(v);
%
%   Input parameters:
%       v   : vector of sample (1 dimention)
%   Output parameters:
%       x   : coordinate along x
%       y   : coordinate along y
%
%   This function compute a discrete continuous density function from a
%   sample list. 
%   
%   Example:
%
%       gsp_point2dcdf([0 3 4 2 2 0 1 0 1 1 2])
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_point2dcdf.php

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

N = numel(v);
tol = 1e-10;
v = sort(v);
[x, inds] = unique(round(v(:)*1/tol)*tol);

y = (inds-1)/(N-1);

%% Old code
%     N = numel(v);
% 
%     v = sort(v(:));
%     
%     x = zeros(N,1);
%     y = zeros(N,1);
%     x(1) = v(1);
%     y(1) = 1;
%     ind = 1;
%     for ii = 1:(N-1)
%         if v(ii+1)>v(ii)+ eps(10);
%             x(ind+1) = v(ii+1);
%             y(ind+1) = y(ind)+1;                
%             ind = ind + 1;
%         else
%             y(ind) = y(ind) + 1;
%             if ii == N-1
%                 x(ind) = v(ii);
%             end
%         end
%         
%     end
%     
%     [~,~,x] = find(x);
%     [~,~,y] = find(y);
%     
%     y = y/N;
%     



end
