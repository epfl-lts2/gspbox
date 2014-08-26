function ldltest
%LDLTEST test program for LDL
%
% Example:
%   ldltest
% See also ldlsparse.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldltest.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

help ldlsparse

A = sparse ([ ]) ;

[L, D, Parent, fl] = ldlsparse (A) ;					    %#ok
[L, D, Parent, fl] = ldlsparse (A, [ ]) ;				    %#ok

try
    [L, D, Parent, fl] = ldlsparse (A, [1 2]) ;
    fprintf ('L = ') ; disp (L)
    fprintf ('D = ') ; disp (D)
    fprintf ('Parent = ') ; disp (Parent)
    fprintf ('fl = %g\n', fl) ;
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok)
    error ('?') ;
end

try
    ldlsparse
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok)
    error ('?')
end

try
    [L, D, Parent, fl] = ldlsparse (1) ;				    %#ok
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok)
    error ('?')
end

try
    x = ldlsparse (1,2) ;						    %#ok
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok)
    error ('?')
end

A =[ ...
1.7     0     0     0     0     0     0     0   .13     0
  0   1.0     0     0   .02     0     0     0     0   .01
  0     0   1.5     0     0     0     0     0     0     0
  0     0     0   1.1     0     0     0     0     0     0
  0   .02     0     0   2.6     0   .16   .09   .52   .53
  0     0     0     0     0   1.2     0     0     0     0
  0     0     0     0   .16     0   1.3     0     0   .56
  0     0     0     0   .09     0     0   1.6   .11     0
.13     0     0     0   .52     0     0   .11   1.4     0
  0   .01     0     0   .53     0   .56     0     0   3.1 ] ;
A = sparse (A) ;
b = [ ...
 .98 .64 .05
 .58 .20 .44
 .42 .37 .30
 .51 .78 .84
 .33 .68 .01
 .43 .46 .76
 .22 .56 .97
 .57 .79 .99
 .76 .05 .78
 .52 .60 .43 ] ;
P = [3 10 2 5 8 6 9 7 1 4] ;
I = speye (10) ;

[L, D, Parent, fl] = ldlsparse (A) ;
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14)
    error ('?') ;
end ;

Parent2 = etree (A) ;
if (any (Parent2 ~= Parent))
    error ('?') ;
end

[L, D, Parent] = ldlsparse (A) ;					    %#ok
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

[L, D] = ldlsparse (A) ;
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

L2 = ldlsparse (A) ;
err = norm (L - L2, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

x = ldlsparse (A, [ ], b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

[x, fl] = ldlsparse (A, [ ], b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14)
    error ('?') ;
end ;

[L, D, Parent, fl] = ldlsparse (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14)
    error ('?') ;
end ;

Parent2 = etree (A (P,P)) ;
if (any (Parent2 ~= Parent))
    error ('?') ;
end

clf
subplot (2,2,1), spy (A),           title ('original matrix') ;
subplot (2,2,2), spy (A (P,P)),     title ('permuted matrix') ;
subplot (2,2,3), spy (L+D+L'),      title ('L+D+L''') ;
subplot (2,2,4), treeplot (Parent), title ('elimination tree') ;

[L, D, Parent] = ldlsparse (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

[L, D] = ldlsparse (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

L2 = ldlsparse (A, P) ;
err = norm (L - L2, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

x = ldlsparse (A, P, b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14)
    error ('?') ;
end ;

[x, fl] = ldlsparse (A, P, b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14)
    error ('?') ;
end ;

fprintf ('\nldl: all tests passed\n') ;

