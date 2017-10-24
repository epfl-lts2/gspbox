% sgwt_cheby_eval : Evaluate shifted Chebyshev polynomial on given domain
%
% function r=sgwt_cheby_eval(x,c,arange)
%
% Compute Chebyshev polynomial of laplacian applied to input.
% This is primarily for visualization
%
% Inputs:
% x - input values to evaluate polynomial on
% c - Chebyshev coefficients (c(1+j) is jth coefficient)
% arange - interval of approximation. Note that x need not be inside
%          arange, but the polynomial will no longer be near the
%          approximated function outside of arange.

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

function r=sgwt_cheby_eval(x,c,arange)
% By setting the operator L to mean (pointwise) multiplication by x,
% and f to be vector of ones, p(L)f will contain p(x) at each
% point. This lets us use sgwt_cheby_op to evaluate the Chebyshev polynomials.
L=spdiags(x(:),0,speye(numel(x)));
f=ones(size(x(:)));
r=sgwt_cheby_op(f,L,c,arange);

if iscell(r)
  for k=1:numel(r)
    r{k}=reshape(r{k},size(x));
  end
else
  r=reshape(r,size(x));
end



  
  
