function C = gsp_strong(A,B)
%GSP_STRONG Strong product between matrices
%   Usage:  C = gsp_strong(A,B)
%
%   Input parameters:
%         A          : matrix size N1*M1
%         B          : matrix size N2*M2
%   Output parameters:
%         C          : matrix size (N1xN2)  (M1xM2)
%
%   'gsp_strong' computes the strong product between two matrices A
%   and B. The cartesian product between two matrices is
%
%      C = strong(A,B) = kron(A,B) + cartesian(A,B) 
%
%   See also: gsp_graph_product
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_strong.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

% Author : Francesco Grassi
% Date   : July 2016

C = kron(A,B) + gsp_cartesian(A,B);


end
