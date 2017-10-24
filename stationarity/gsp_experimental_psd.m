function gf = gsp_experimental_psd(G,C)
%GSP_EXPERIMENTAL_PSD Experimental power density function
%   Usage:  f_hat=gsp_gft(G,f);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         C          : Covariance matrix
%   Output parameters:
%         gf         : PSD filter
%
%   This function estimate the PSD from the covariance matrix with
%
%      T = U' * C * U 
%
%   where U is the graph Fourier basis. The function then interpolates
%   the diagonal of T with splines.
%
%   To compute the Fourier basis of a graph G, you can use the function:
%
%           G = gsp_compute_fourier_basis(G);
%
%
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_experimental_psd.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781


% Author : Nathanael Perraudin
% Date: 6 January 2016


if ~isnumeric(G);
    if ~isfield(G,'U')
       error(['GSP_EXPERIMENTAL_PSD: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    U = G.U;
else
    U = G;
end


    Cf = U'*C*U;
    
    cf = real(diag(Cf));
    [ e,ind ]= unique(G.e);
    gf = @(s) max(spline(e,cf(ind),s),0);
end

% % Old code with polynomial fitting
%     X = zeros(N,Order+1);
%     X(1,:) = 1;
%     for ii = 1:Order
%         X(ii+1) = X(ii+1) .* G.e;
%     end
%     
%     alpha = Pinv(X)*cf;
%     [p,S,mu] = polyfit(G.e,cf,Order);
%     gf = @(x) abs(polyval(p,x,S,mu))+ eps;


