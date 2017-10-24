% sgwt_soft_threshold : Soft thresholding operator
%
% x_t = bpdq_soft_threshold(x,tgamma)
%
% Applies soft thresholding  to each component of x
%
% Inputs:
% x - input signal
% tgamma - threshold
%
% Outputs:
% x_t - soft thresholded result

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.


function x_t = sgwt_soft_threshold(x,tgamma)
  tmp=abs(x)-tgamma;
  x_t = sign(x).*tmp.*(tmp>0);
  
