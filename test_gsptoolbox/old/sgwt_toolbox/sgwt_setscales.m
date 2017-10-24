% sgwt_setscales : Compute a set of wavelet scales adapted to spectrum bounds
%
% function s=sgwt_setscales(lmin,lmax,Nscales)
% 
% returns a (possibly good) set of wavelet scales given minimum nonzero and 
% maximum eigenvalues of laplacian
% 
% returns scales logarithmicaly spaced between minimum and maximum
% "effective" scales : i.e. scales below minumum or above maximum
% will yield the same shape wavelet (due to homogoneity of sgwt kernel : 
% currently assuming sgwt kernel g given as abspline with t1=1, t2=2)
%
% Inputs : 
% lmin,lmax - minimum nonzero and maximum eigenvalue of
%             Laplacian. Note that in design of transform with
%             scaling function, lmin may be taken just as a fixed
%             fraction of lmax,  and may not actually be the
%             smallest nonzero eigenvalue 
% Nscales - # of wavelet scales
%
% Outputs :
% s - wavelet scales

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


function s=sgwt_setscales(lmin,lmax,Nscales)
  t1=1;
  t2=2;
  
  smin=t1/lmax;
  smax=t2/lmin;
  % scales should be decreasing ... higher j should give larger s
  s=exp(linspace(log(smax),log(smin),Nscales));
  
