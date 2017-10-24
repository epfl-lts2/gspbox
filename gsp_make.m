function [ ] = gsp_make( )
%GSP_MAKE Compile the necessary toolboxes for the gspbox
%   Usage: gsp_make();
%
%   This function compile the routine for the gspbox::
%       
%           gsp_make();
%

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
