% sgwt_framebounds : Compute approximate frame bounds for given sgw transform
%
% function [A,B,sg2,x]=sgwt_framebounds(g,lmin,lmax)
%
% Inputs : 
% g - function handles computing sgwt scaling function and wavelet
% kernels
% lmin,lmax - minimum nonzero, maximum eigenvalue
%
% Outputs :
% A , B - frame bounds
% sg2 - array containing sum of g(s_i*x)^2 (for visualization)
% x - x values corresponding to sg2

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

function [A,B,sg2,x]=sgwt_framebounds(g,lmin,lmax)
  N=1e4; % number of points for line search
  x=linspace(lmin,lmax,N);
  Nscales=numel(g);

  sg2=zeros(size(x));
  for ks=1:Nscales
    sg2=sg2+(g{ks}(x)).^2;
  end
  A=min(sg2);
  B=max(sg2);
  
  
