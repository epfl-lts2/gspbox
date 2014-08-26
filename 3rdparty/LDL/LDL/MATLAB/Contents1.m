% LDL package: simple sparse LDL factorization
%
% Primary routines:
% 
%   ldlsparse   - LDL' factorization of a real, sparse, symmetric matrix
%   ldlsymbol   - symbolic Cholesky factorization
%
% Helper routines:
%
%   ldldemo     - demo program for LDL
%   ldlrow      - an m-file description of the algorithm used by LDL
%   ldltest     - test program for LDL
%   ldlmain2    - compiles and runs a longer test program for LDL
%   ldl_install - compile and install the LDL package for use in MATLAB.
%   ldl_make    - compile LDL
%
% Example:
%
%       [L, D, Parent, fl] = ldlsparse (A)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/Contents1.php

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

% LDL License:  GNU Lesser General Public License as published by the
%   Free Software Foundation; either version 2.1 of the License, or
%   (at your option) any later version.
%
% Acknowledgements:
%
%   This work was supported by the National Science Foundation, under
%   grant CCR-0203270.
%
%   Portions of this work were done while on sabbatical at Stanford University
%   and Lawrence Berkeley National Laboratory (with funding from the SciDAC
%   program).  I would like to thank Gene Golub, Esmond Ng, and Horst Simon
%   for making this sabbatical possible.

