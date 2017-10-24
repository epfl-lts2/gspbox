function [  ] = gsp_install_unlocbox(  )
%GSP_INSTALL_UNLOCBOX Install the UNLocBoX from the web
%   Usage: gsp_install_unlocbox();
%  
%   This function installs the UNLocBoX. It require an internet connection.
%
global GLOBAL_gsppath;

FS=filesep;

fprintf('Install the UNLocBoX: ')
try
    outputdir = [GLOBAL_gsppath,FS,'3rdparty'];
    if isunix
        tarfilename = webread('https://github.com/epfl-lts2/unlocbox/releases/download/1.7.4/unlocbox-1.7.4.tar.gz');
        untar(tarfilename,outputdir)        
    else
        zipfilename = webread('https://github.com/epfl-lts2/unlocbox/releases/download/1.7.4/unlocbox-1.7.4.zip');
        unzip(zipfilename,outputdir)
    end
    fprintf('Installation sucessfull!\n\n')
catch
    warning('Could not install the UNLocBox')
end


end
