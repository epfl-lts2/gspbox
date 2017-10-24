% sgwt_forward : Compute forward of sgw transform
%
% function y=sgwt_forward(f,L,c,arange)
%
% This routine is just a wrapper calling sgwt_cheby_op, and is included for 
% convenience
% 
% Inputs:
% f - Input signal
% L - laplacian
% c - cell array of Chebyshev coefficients defining transform
% arange - spectral approximation range
%
% Outputs:
% y - SGWT coefficients

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
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.function r=sgwt_forward(f,L,c,arange)
r=sgwt_cheby_op(f,L,c,arange);
