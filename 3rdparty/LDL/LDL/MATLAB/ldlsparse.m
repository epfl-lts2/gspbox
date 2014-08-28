function [arg1, arg2, arg3, arg4] = ldlsparse (A, P, b)			    %#ok
%LDLSPARSE LDL' factorization of a real, sparse, symmetric matrix
%
% Example:
%       [L, D, Parent, fl] = ldlsparse (A)
%       [L, D, Parent, fl] = ldlsparse (A, P)
%       [x, fl] = ldlsparse (A, [ ], b)
%       [x, fl] = ldlsparse (A, P, b)
%
% Let I = speye (size (A,1)). The factorization is (L+I)*D*(L+I)' = A or A(P,P).
% A must be sparse, square, and real.  Only the diagonal and upper triangular
% part of A or A(P,P) are accessed.  L is lower triangular with unit diagonal,
% but the diagonal is not returned.  D is a diagonal sparse matrix.  P is either
% a permutation of 1:n, or an empty vector, where n = size (A,1).  If not
% present, or empty, then P=1:n is assumed.  Parent is the elimination tree of
% A or A(P,P).  If positive, fl is the floating point operation count, or
% negative if any entry on the diagonal of D is zero.
%
% In the x = ldlsparse (A, P, b) usage, the LDL' factorization is not returned.
% Instead, the system A*x=b is solved for x, where both b and x are dense.
%
% If a zero entry on the diagonal of D is encountered, the LDL' factorization is
% terminated at that point.  If there is no fl output argument, an error occurs.
% Otherwise, fl is negative, and let d=-fl.  D(d,d) is the first zero entry on
% the diagonal of D.  A partial factorization is returned.  Let B = A, or A(P,P)
% if P is present.  Let F = (L+I)*D*(L+I)'.  Then F (1:d,1:d) = B (1:d,1:d).
% Rows d+1 to n of L and D are all zero.
%
% See also chol, ldl, ldlsymbol, symbfact, etree
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldlsparse.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

error ('ldlsparse mexFunction not found') ;

