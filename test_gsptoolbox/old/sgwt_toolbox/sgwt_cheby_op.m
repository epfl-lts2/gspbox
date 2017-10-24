% sgwt_cheby_op : Chebyshev polynomial of Laplacian applied to vector
%
% function r=sgwt_cheby_op(f,L,c,arange)
%
% Compute (possibly multiple) polynomials of laplacian (in Chebyshev
% basis) applied to input.
%
% Coefficients for multiple polynomials may be passed as a cell array. This is
% equivalent to setting
% r{1}=sgwt_cheby_op(f,L,c{1},arange);
% r{2}=sgwt_cheby_op(f,L,c{2},arange);
% ...
% 
% but is more efficient as the Chebyshev polynomials of L applied
% to f can be computed once and shared.
%
% Inputs:
% f- input vector
% L - graph laplacian (should be sparse)
% c - Chebyshev coefficients. If c is a plain array, then they are
%     coefficients for a single polynomial. If c is a cell array, 
%     then it contains coefficients for multiple polynomials, such 
%     that c{j}(1+k) is k'th Chebyshev coefficient the j'th polynomial.
% arange - interval of approximation
%
% Outputs:
% r - result. If c is cell array, r will be cell array of vectors
%     size of f. If c is a plain array, r will be a vector the size
%     of f.

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

function r=sgwt_cheby_op(f,L,c,arange)

if ~iscell(c)
  r=sgwt_cheby_op(f,L,{c},arange);
  r=r{1};
  return;
end

Nscales=numel(c);
M=zeros(size(Nscales));
for j=1:Nscales
    M(j)=numel(c{j});
end
assert(all(M>=2));

maxM=max(M);
%Twf_new = T_j(L) f
%Twf_cur T_{j-1}(L) f
%TWf_old T_{j-2}(L) f

a1=(arange(2)-arange(1))/2;
a2=(arange(2)+arange(1))/2;

Twf_old=f; %j=0;
Twf_cur=(L*f-a2*f)/a1; % j=1;
for j=1:Nscales
    r{j}=.5*c{j}(1)*Twf_old + c{j}(2)*Twf_cur;
end

for k=2:maxM
    Twf_new = (2/a1)*(L*Twf_cur-a2*Twf_cur)-Twf_old;
    for j=1:Nscales
        if 1+k<=M(j)
            r{j}=r{j}+c{j}(k+1)*Twf_new;
        end
    end
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
end
