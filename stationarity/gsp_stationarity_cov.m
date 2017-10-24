function C = gsp_stationarity_cov(X)
%GSP_STATIONARITY_COV Covariance matrix from graph stationary data
%   Usage:  C = gsp_stationarity_cov(X)
%
%   Input parameters:
%         X          : Data (M x N matrix) 
%   Output parameters:
%         C          : Covariance matrix (M x M)
%
%   This function estimates the covariance from the data. Every sample has
%   the same expected average.
%
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_stationarity_cov.html

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

Nso = size(X,2);
X = X - mean(X(:));
C = (X*X')/Nso;


end
