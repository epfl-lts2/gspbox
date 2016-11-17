function gsp_start()
%GSP_START Initialize the toolbox
%   Usage: gsp_start();
%
%   Initialisation script for the GSPBox. This script add the different
%   path needed to run the toolbox. 
%
%   References:
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/gsp_start.php

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


% Author: Nathanael Perraudin
% Date: 14 March 2014


%% adding dependency

global GLOBAL_gsppath;
GLOBAL_gsppath = fileparts(mfilename('fullpath'));


addpath(genpath(GLOBAL_gsppath));

% Load the version number
bp=[GLOBAL_gsppath,filesep];
[FID, MSG] = fopen ([bp,'gspbox_version'],'r');
if FID == -1
    error(MSG);
else
    gspbox_version = fgetl (FID);
    fclose(FID);
end

banner = sprintf(['GSPBox version %s. Copyright 2013-2015 LTS2-EPFL,\n',...
                  'by Nathanael Perraudin, Johan Paratte, David Shuman ',...
                  'and Vassilis Kalofolias'], ...
                   gspbox_version);

% display banner
disp(banner);


end

