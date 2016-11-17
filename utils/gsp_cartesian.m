function C = gsp_cartesian(A,B)
%GSP_CARTESIAN Cartesian product between vectors or matrices
%   Usage:  C = gsp_cartesian(A,B)
%
%   Input parameters:
%         A          : vector size N1 or matrix size N1 x M1
%         B          : vector size N2 or matrix size N2 x M2
%   Output parameters:
%         C          : vector size N1*N2 or matrix size (N1*N2) x (M1*M2)
%
%   If A and B are vectors 'gsp_cartesian' computes the following product:
%
%      C = cartesian(A,B) = kron(A,ones(size(B))) + kron(ones(size(A),B)
%
%   If A and B are matrices 'gsp_cartesian' computes the following product:
%
%      C = cartesian(A,B) = kron(A,eye(size(B))) + kron(eye(size(A),B)
%
%   See also: gsp_graph_product 
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_cartesian.php

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

%TO DO: extend to more than 2 factors

if isvector(A) && isvector(B)
    C = kron(A,ones(size(B))) + kron(ones(size(A)),B);
elseif ismatrix(A) && ismatrix(B)
    C = kron(A,eye(size(B))) + kron(eye(size(A)),B);
else
    error('GSP_CARTESIAN: Cartesian product can be computed only between matrices or between vectors')
end

end
