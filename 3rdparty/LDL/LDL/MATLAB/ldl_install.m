function ldl_install
%LDL_INSTALL compile and install the LDL package for use in MATLAB.
% Your current working directory must be LDL for this function to work.
%
% Example:
%   ldl_install
%
% See also ldlsparse, ldlsymbol
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldl_install.php

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

ldl_make
addpath (pwd) ;
fprintf ('LDL has been compiled and installed.  The path:\n') ;
disp (pwd) ;
fprintf ('has been added to your path.  Use pathtool to add it permanently.\n');
ldldemo

