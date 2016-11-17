function [y,x] = gsp_ddf2dcdf(v,x)
%GSP_DDF2DCDF Discrete density function to dicrete cumulative density function
%   Usage : [y] = gsp_ddf2dcdf(v);
%           [y,x] = gsp_ddf2dcdf(v,x);
%
%   Input parameters:
%       v   : Discrete density function
%       x   : Axis value (Optional)
%   Output parameters:
%       y   : Discrete cumulative ensity function
%       x   : Axis value (Same as input)
%
%   This function goes from discrete density function to discrete
%   cumulative density function.
%   
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_ddf2dcdf.php

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


N = length(v);
y = zeros(N,1);
y(1) = v(1);
for ii = 1:(N-1) 
    y(ii+1) = y(ii) +v(ii+1);
end

y = y/y(end);

end
