function y = gsp_smooth_step(x, a)
%GSP_SMOOTH_STEP Smooth step function from 0 to 1 arround 0.5
%   Usage: y = gsp_smooth_step(x);
%          y = gsp_smooth_step(x, a)
%   
%   Input parameters:
%          x        : input value
%          a        : smoothing parameter (default 1)
%   Ouput parameters:
%          y        : output value of the function
%
%   This function is a smooth step from 0 to 1 arround 0 using the
%   following function:
%
%               /   0                                      if x < 0
%       s(x) = | exp(-a/x) / ( exp(-a/x) + exp(-a/(1-x)) ) if x in [0, 1]
%               \   1                                      if x > 1
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_smooth_step.php

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
% Date  : 5 mai 2016


if nargin<2
    a = 1;
end

y=fx(x,a);
y = y./(y+fx(1-x,a));

end

function y=fx(x,a)
    y = exp(-a./x);
    y(x<0)=0;
end

