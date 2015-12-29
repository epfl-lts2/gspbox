function [L, D] = ldlrow (A)
%LDLROW an m-file description of the algorithm used by LDL
%
% Example:
%  [L, D] = ldlrow (A)
%
%  Compute the L*D*L' factorization of A, by rows.  Returns
%  full L and D matrices.  This routine serves as an outline
%  of the numerical factorization performed by ldl.c.
%
%  Here is a diagram of how L is computed.  "a" means an
%  entry that is accessed at the kth step, and "c" means an
%  entry that is computed.  A "-" means neither accessed nor
%  computed.  A "1" means the value of the entry is L (the
%  unit diagonal of L), and it is accessed at the kth step.
%  A "." means the value is zero.
%
%  The L matrix
%
%     1 . . . . . . .
%     a 1 . . . . . .
%     a a 1 . . . . .
%     a a a 1 . . . .
%     c c c c c . . .  <- kth row of L
%     - - - - - - . .
%     - - - - - - - .
%     - - - - - - - -
%
%  The A matrix:
%
%             the kth column of A
%             v
%     - - - - a - - -
%     - - - - a - - -
%     - - - - a - - -
%     - - - - a - - -
%     - - - - a - - -  <- kth row of A
%     - - - - - - - - 
%     - - - - - - - -
%     - - - - - - - -
%
%  The D matrix:
%
%             the kth column of D
%             v
%     a . . . . . . .
%     . a . . . . . .
%     . . a . . . . .
%     . . . a . . . .
%     . . . . c . . .  <- kth row of D
%     . . . . . . . . 
%     . . . . . . . .
%     . . . . . . . .
%
% See also ldlsparse.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldlrow.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Copyright 2006-2007 by Timothy A. Davis, http://www.suitesparse.com

[m n] = size (A) ;
L = zeros (n, n) ;
D = zeros (n, 1) ;
A = full (A) ;

L (1, 1) = 1 ;
D (1) = A (1,1) ;

for k = 2:n

    % note the sparse triangular solve.  For the sparse
    % case, the pattern of y is the same as the pattern of
    % the kth row of L.
    y = L (1:k-1, 1:k-1) \ A (1:k-1, k) ;

    % scale row k of L
    L (k, 1:k-1) = (y ./ D (1:k-1))' ;
    L (k, k) = 1 ;

    % compute the diagonal
    D (k) = A (k,k) - L (k, 1:k-1) * y ;
end

D = diag (D) ;

