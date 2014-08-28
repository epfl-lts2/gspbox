function [  ] = gsp_install(  )
%GSP_INSTALL Install third party software
%   Usage: gsp_install();
%  
%   This function install third party software. It require an internet
%   connection.
%
%   It will install the gaimc toolbox, the UNLocBoX and compile some
%   functions.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/gsp_install.php

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


global GLOBAL_gsppath;
gsp_start;

FS=filesep;

fprintf('Install the GAIMC toolbox: ')
try
    zipfilename = 'http://www.mathworks.com/matlabcentral/fileexchange/downloads/114153/download';
    outputdir = [GLOBAL_gsppath,FS,'3rdparty',FS,'gaimc'];
    unzip(zipfilename,outputdir)
    fprintf('Installation sucessfull!\n\n')
catch
    warning('Could not install the UNLocBoX')
end

fprintf('Install the UNLocBoX: ')
try
    outputdir = [GLOBAL_gsppath,FS,'3rdparty'];
    if isunix
        tarfilename = 'http://unlocbox.sourceforge.net/unlocbox.tar.gz';
        untar(tarfilename,outputdir)        
    else
        zipfilename = 'http://unlocbox.sourceforge.net/unlocbox.zip';
        unzip(zipfilename,outputdir)
    end
    fprintf('Installation sucessfull!\n\n')
catch
    warning('Could not install the UNLocBox')
end

gsp_make;


end


