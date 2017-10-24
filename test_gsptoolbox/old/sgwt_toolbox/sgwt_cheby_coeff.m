% sgwt_cheby_coeff : Compute Chebyshev coefficients of given function
%
% function c=sgwt_cheby_coeff(g,m,N,arange)
%
% Inputs:
% g - function handle, should define function on arange
% m - maximum order Chebyshev coefficient to compute
% N - grid order used to compute quadrature (default is m+1)
% arange - interval of approximation (defaults to [-1,1] )
% 
% Outputs:
% c - array of Chebyshev coefficients, ordered such that c(j+1) is 
% j'th Chebyshev coefficient

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

function c=sgwt_cheby_coeff(g,m,N,arange)
  if ~exist('N','var')
    N=m+1;
  end
  if ~exist('arange','var')
      arange=[-1, 1];
  end
  a1=(arange(2)-arange(1))/2;
  a2=(arange(2)+arange(1))/2;
  for j=1:m+1
    c(j)=sum (g(a1* cos( (pi*((1:N)-0.5))/N) + a2).*cos(pi*(j-1)*((1:N)-.5)/N) )*2/N;
  end
  
