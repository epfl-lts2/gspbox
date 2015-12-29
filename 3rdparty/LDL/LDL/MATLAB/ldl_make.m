function ldl_make
%LDL_MAKE compile LDL
%
% Example:
%   ldl_make       % compiles ldlsparse and ldlsymbol
%
% See also ldlsparse, ldlsymbol
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldl_make.php

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

d = '' ;
if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims' ;
end
eval (sprintf ('mex -O %s -DLDL_LONG -I../../SuiteSparse_config -I../Include -output ldlsparse ../Source/ldl.c ldlmex.c', d)) ;
eval (sprintf ('mex -O %s -DLDL_LONG -I../../SuiteSparse_config -I../Include -output ldlsymbol ../Source/ldl.c ldlsymbolmex.c', d)) ;

eval (sprintf ('mex -O %s -DLDL_LONG -I../../SuiteSparse_config -I../Include -output ldlsparse_short ../Source/ldl.c ldlmex_short.c', d)) ;
eval (sprintf ('mex -O %s -DLDL_LONG -I../../SuiteSparse_config -I../Include -output ldlsymbol_extra ../Source/ldl.c ldlsymbol_extra.c', d)) ;
eval (sprintf ('mex -O %s -DLDL_LONG -I../../SuiteSparse_config -I../Include -output ldlnumeric ../Source/ldl.c ldlnumericmex.c', d)) ;
fprintf ('LDL successfully compiled.\n') ;


