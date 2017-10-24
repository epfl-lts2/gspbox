% sgwt_inverse : Compute inverse sgw transform, via conjugate gradients
%
% function r=sgwt_inverse(y,L,c,arange)
%
% Inputs:
% y - sgwt coefficients
% L - laplacian
% c - cell array of Chebyshev coefficients defining transform
% arange - spectral approximation range
%
% Selectable Control Parameters
% tol - tolerance for conjugate gradients (default 1e-6)
%

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

function r=sgwt_inverse(y,L,c,arange,varargin)
control_params={'tol',1e-6};
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

assert(iscell(c));
N=size(L,1);
% first compute adj = W^*y ( sort of slowly )
fprintf('computing adjoint\n');
adj=sgwt_adjoint(y,L,c,arange);
% W^* W
% compute P(x) = p(x)^2
fprintf('computing cheby coeff for P=p^2\n');
for j=1:numel(c)
    M(j)=numel(c{j});
end
maxM=max(M);
% dkh : code below could remove unnecessary use of cell arrays.
d{1}=zeros(1,1+2*(maxM-1));
for j=1:numel(c)
    cpad{j}=zeros(maxM,1);
    cpad{j}(1:M(j))=c{j};
    d{1}=d{1}+sgwt_cheby_square(cpad{j});
end
wstarw = @(x) sgwt_cheby_op(x,L,d{1},arange);
%% conjugate gradients
fprintf('computing inverse by conjugate gradients\n');
r=pcg(wstarw,adj,tol);
