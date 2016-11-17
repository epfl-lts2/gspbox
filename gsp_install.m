function [  ] = gsp_install( install_unlcobox  )
%GSP_INSTALL Install third party software
%   Usage: gsp_install();
%  
%   Input parameters
%       install_unlcobox    : bool (1 to install the UNLocBoX)
%
%   This function installs third party software. It require an internet
%   connection.
%
%   You should install the UNLocBoX only if you are not using the
%   developement version.
%
%   It will install the gaimc toolbox and compile some functions.
%
%   See also: gsp_install_unlocbox
%
%   Url: http://lts2research.epfl.ch/gsp/doc/gsp_install.php

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


global GLOBAL_gsppath;
gsp_start;

if nargin<1
    install_unlcobox = 0;
end

FS=filesep;

fprintf('Install the GAIMC toolbox: ')
try
    zipfilename = 'http://www.mathworks.com/matlabcentral/fileexchange/submissions/24134/v/1/download/zip';
    outputdir = [GLOBAL_gsppath,FS,'3rdparty',FS,'gaimc'];
    unzip(zipfilename,outputdir)
    fprintf('Installation sucessfull!\n\n')
catch
    warning('Could not install GAIMC')
end

if install_unlcobox
    gsp_install_unlocbox
else
    fprintf('The function did not install the UNLocBoX. To do it, run: gsp_install_unlcobox')
end

gsp_make;


end


