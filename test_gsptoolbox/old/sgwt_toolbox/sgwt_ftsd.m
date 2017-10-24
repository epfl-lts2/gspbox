% sgwt_ftsd : Compute forward transform in spectral domain
% 
% function r=sgwt_ftsd(f,g,t,L)
%
% Compute forward transform by explicitly computing eigenvectors and 
% eigenvalues of graph laplacian
%
% Uses persistent variables to store eigenvectors, so decomposition
% will be computed only on first call
%
% Inputs:
% f - input data
% g - sgw kernel
% t - desired wavelet scale
% L - graph laplacian
%
% Outputs:
% r - output wavelet coefficients

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

function r=sgwt_ftsd(f,g,t,L)
persistent V D Lold
if (isempty(V) || any(vec(L~=Lold)))
  fprintf('Diagonalizing %g x %g L (could take some time ...)\n',size(L,1),size(L,2));
  [V,D]=eig(full(L));
  Lold=L;
end
lambda=diag(D);
fhat=V'*f;
r=V*(fhat.*g(t*lambda));
  
