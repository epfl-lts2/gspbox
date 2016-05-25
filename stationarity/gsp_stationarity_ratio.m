function r = gsp_stationarity_ratio(G, C)
%GSP_STATIONARITY_RATIO Assert the stationarity level of some data
%   Usage:  r = gsp_stationarity_ratio(G, C)
%
%   Input parameters:
%         G          : Graph
%         C          : Covariance matrix
%   Output parameters:
%         r          : Ratio
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
%   Url: http://lts2research.epfl.ch/gsp/doc/stationarity/gsp_stationarity_ratio.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

CF = G.U' * C * G.U;

r = gsp_diagonal_ratio(CF);



end


function r = gsp_diagonal_ratio(M)

if size(M,1) == 1 || size(M,2) ==1
    error('This function acts on martrices.')
end

r = norm(diag(M))/norm(M,'fro');

end
