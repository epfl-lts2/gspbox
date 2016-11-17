function [ ] = gsp_make( )
%GSP_MAKE Compile the necessary toolboxes for the gspbox
%   Usage: gsp_make();
%
%   This function compile the routine for the gspbox:
%       
%           gsp_make();
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/gsp_make.php

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

% TODO: clean this function!

gsp_start;
global GLOBAL_gsppath;

FS=filesep;

test = 0 ;
paths = { } ;
gsp_path = GLOBAL_gsppath ;
paths = add_to_path (paths, gsp_path) ;

% compile and install AMD
try
    paths = add_to_path (paths, [gsp_path,FS,'3rdparty',FS,'LDL',FS,'AMD',FS,'MATLAB']) ;
    amd_make ;
catch me
    disp (me.message) ;
    fprintf ('AMD not installed\n') ;
end

% compile and install LDL
try
    paths = add_to_path (paths, [gsp_path,FS,'3rdparty',FS,'LDL',FS,'LDL',FS,'MATLAB']) ;
    ldl_make ;
    if test
        ldlmain2 ;
        ldltest ;  
    end
catch me
    disp (me.message) ;
    fprintf ('LDL not installed\n') ;
end

cd (gsp_path)


end

function paths = add_to_path (paths, newpath)
    % add a path
    cd (newpath) ;
    addpath (newpath) ;
    paths = [paths { newpath } ] ;						 
end

