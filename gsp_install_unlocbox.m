function [  ] = gsp_install_unlocbox(  )
%GSP_INSTALL_UNLOCBOX Install the UNLocBoX from the web
%   Usage: gsp_install_unlocbox();
%  
%   This function installs the UNLocBoX. It require an internet connection.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/gsp_install_unlocbox.php

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

FS=filesep;

fprintf('Install the UNLocBoX: ')
try
    outputdir = [GLOBAL_gsppath,FS,'3rdparty'];
    if isunix
        tarfilename = webread('http://lts2.epfl.ch/unlocbox/download_unix.php');
        untar(tarfilename,outputdir)        
    else
        zipfilename = webread('https://lts2.epfl.ch/unlocbox/download_win.php');
        unzip(zipfilename,outputdir)
    end
    fprintf('Installation sucessfull!\n\n')
catch
    warning('Could not install the UNLocBox')
end


end

