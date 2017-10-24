% sgwt_adjoint : Compute adjoint of sgw transform
%
% function adj=sgwt_inverse(y,L,c,arange)
%
% Inputs:
% y - sgwt coefficients
% L - laplacian
% c - cell array of Chebyshev coefficients defining transform
% arange - spectral approximation range
%
% Outputs:
% adj - computed sgwt adjoint applied to y

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

function adj=sgwt_adjoint(y,L,c,arange)

assert(iscell(c));
N=size(L,1);
% first compute adj = W^*y ( sort of slowly )
adj=zeros(N,1);
%fprintf('computing adjoint\n');
for j=1:numel(c)
    tmp=sgwt_cheby_op(y{j},L,c{j},arange);
    adj=adj+tmp;
end

