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

