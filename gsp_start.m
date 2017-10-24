function gsp_start()
%GSP_START Initialize the toolbox
%   Usage: gsp_start();
%
%   Initialisation script for the GSPBox. This script add the different
%   path needed to run the toolbox. 
%
%   References: perraudin2014gspbox
%


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
