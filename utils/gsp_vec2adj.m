function W = gsp_vec2adj(A,x,type)
%GSP_VEC2ADJ Create the matrix W with sparsity pattern A and entries x
%
%   Usage: W = gsp_vec2adj(A,x)
%              gsp_vec2adj(A,x,type)
%
%   Input parameters:
%       A       : Matrix of sparsity pattern (e.g. binary adjacency matrix)
%       x       : Vector of entries
%       type    : Type of matrix (string) (default 'sym' if A is symmetric)
%                      sym: W will be symmetric (x must be of size nnz(A)/2)
%                     asym: W is asymmetric (x must be of size nnz(A))
%   Output parameters
%       W       : Weighted matrix
%
%   Create the matrix W with sparsity pattern A and entries x%
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_vec2adj.php

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

% Author: Francesco Grassi
% Date   : July 2016


if nargin<3
    if issymmetric(full(double(A)))
        type='sym';
    else
        type='asym';
    end
end

if ~nnz(A)==length(x) && ~nnz(A)/2==length(x)
    error('x must have the same size of nnz(A) or nnz(A)/2 if A is symmetric')
end

[N1,N2] = size(A);

switch type
    case 'sym'
        [i,j] = find(triu(A));
        W = sparse(i,j,x,N1,N2)+sparse(i,j,x,N1,N2)';
    case 'asym'
        [i,j] = find(A);
        W = sparse(i,j,x,N1,N2);
    otherwise
        error('Unknow type');
end

    




