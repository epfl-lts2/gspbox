function r = gsp_stationarity_ratio(G, C, param)
%GSP_STATIONARITY_RATIO This is a number from [0 to 1] depicting how close to the PCA basis the graph basis is
%
%   The method examines the percentage of the data variance that is not in
%   the diagonal of the covariance of the GFT of X. The index can be used to
%   describe how well a graph fits a given data matrix X, or distribution. 
%   An index of 0 means that the data are graph stationary on G.
%
%   Usage:  r = gsp_stationarity_ratio(G, C)
%
%   Input parameters:
%         G          : Graph
%         C          : Covariance matrix
%   Output parameters:
%         r          : Ratio
%
%   Optional parameters: 
%   params.verbose   :0 = nothing, 1 = plot the covariance of GFT(x) (default 0)
% 
%   This function compute the ratio of energy contained into the diagonal
%   of the Fourier covariance matrix:
%
%      T = U' * C * U 
%
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_stationarity_ratio.html

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


% Authors : Nathanael Perraudin, Andreas Loukas
% Date    : 6 January 2016

if nargin < 3, param = struct(); end

if ~isfield(param,'verbose'), param.verbose = 0; end

if not(isfield(G, 'U'))
    G = gsp_compute_fourier_basis(G);
end


CF = G.U' * C * G.U;

r = gsp_diagonal_ratio(CF);


if param.verbose > 0
    figure; imagesc(abs(CF)); 
    title('G.U^*  * (X * X^*) * G.U')
end

end


function r = gsp_diagonal_ratio(M)

if size(M,1) == 1 || size(M,2) ==1
    error('This function acts on martrices.')
end

r = norm(diag(M))/norm(M,'fro');

end
