function ldlmain2
%LDLMAIN2 compiles and runs a longer test program for LDL
%
% Example:
%   ldlmain2
%
% See also ldlsparse.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/3rdparty/LDL/LDL/MATLAB/ldlmain2.php

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

help ldlmain2

detail = 0 ;

ldl_make

d = '' ;
if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims' ;
end

mx = sprintf (...
    'mex -O %s -DLDL_LONG -DDLONG -I../../SuiteSparse_config -I../Include', d) ;

% compile ldlmain without AMD
cmd = sprintf ('%s ../Demo/ldlmain.c ../Source/ldl.c', mx) ;
if (detail)
    fprintf ('%s\n', cmd) ;
end
eval (cmd) ;
ldlmain

% compile ldlamd (ldlmain with AMD)
cmd = sprintf ('%s -I../../AMD/Include', mx) ;

files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (files)
    cmd = sprintf ('%s ../../AMD/Source/%s.c', cmd, files {i}) ;
end

cmd = [cmd ' -DUSE_AMD -output ldlamd ../Demo/ldlmain.c ../Source/ldl.c'] ;
if (detail)
    fprintf ('%s\n', cmd) ;
end
eval (cmd) ;
ldlamd


