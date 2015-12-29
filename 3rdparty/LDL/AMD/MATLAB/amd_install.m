function amd_install
%AMD_INSTALL compile and install amd2 for use in MATLAB
%   Your current directory must be AMD/MATLAB for this function to work.
%
% Example:
%   amd_install
%
% See also amd, amd2.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/AMD/MATLAB/amd_install.php

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

% Copyright 1994-2007, Tim Davis, Patrick R. Amestoy, and Iain S. Duff. 

% This orders the same matrix as the ANSI C demo, amd_demo.c.  It includes an

amd_make
addpath (pwd)
fprintf ('\nThe following path has been added.  You may wish to add it\n') ;
fprintf ('permanently, using the MATLAB pathtool command.\n') ;
fprintf ('%s\n\n', pwd) ;
amd_demo

