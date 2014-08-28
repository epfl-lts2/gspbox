function [Lnz, Parent, fl] = ldlsymbol (A, P)				    %#ok
%LDLSYMBOL symbolic Cholesky factorization
%
% Example:
%       [Lnz, Parent, fl] = ldlsymbol (A)
%       [Lnz, Parent, fl] = ldlsymbol (A, P)
%
% P is a permutation of 1:n, an output of AMD, SYMAMD, or SYMRCM, for example.
% Only the diagonal and upper triangular part of A or A(P,P) is accessed; the
% lower triangular part is ignored.  If P is not provided, then P = 1:n is
% assumed.
%
% The elimination tree is returned in the Parent array.  The number of nonzeros
% in each column of L is returned in Lnz.  fl is the floating point operation
% count for a subsequent LDL' factorization.  This mexFunction replicates the
% following MATLAB computations, using ldl_symbolic:
%
%       Lnz = symbfact (A) - 1 ;
%       Parent = etree (A) ;
%       fl = sum (Lnz . (Lnz + 2)) ;
%
% or, if P is provided,
%
%       Lnz = symbfact (A (P,P)) - 1 ;
%       Parent = etree (A (P,P)) ;
%       fl = sum (Lnz . (Lnz + 2)) ;
%
% Note that this routine is not required by LDL, since LDL does its own
% symbolic factorization.
%
% See also ldlsparse, symbfact, etree
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldlsymbol.php

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

error ('ldlsymbol mexFunction not found') ;

