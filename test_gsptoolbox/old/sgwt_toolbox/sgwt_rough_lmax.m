% sgwt_rough_lmax : Rough upper bound on maximum eigenvalue of L
% 
% function lmax=sgwt_rough_lmax(L)
% 
% Runs Arnoldi algorithm with a large tolerance, then increases
% calculated maximum eigenvalue by 1 percent. For much of the SGWT
% machinery, we need to approximate the wavelet kernels on an
% interval that contains the spectrum of L. The only cost of using
% a larger interval is that the polynomial approximation over the
% larger interval may be a slightly worse approxmation on the
% actual spectrum. As this is a very mild effect, it is not likely
% necessary to obtain very tight bonds on the spectrum of L
%
% Inputs : 
% L - input graph Laplacian
%
% Outputs :
% lmax - estimated upper bound on maximum eigenvalue of L

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

function lmax=sgwt_rough_lmax(L)
opts=struct('tol',5e-3,'p',10,'disp',0);
lmax=eigs(L,1,'lm',opts);
lmax=lmax*1.01; % just increase by 1 percent to be robust to error
