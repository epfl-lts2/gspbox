% sgwt_compile_mex : set environment variables and call makefile to compile 
tmp_dir=pwd;
cd([SGWT_ROOT,'mex']);
setenv('matlabroot',matlabroot);
setenv('mexext',mexext);
!make clean
!make
cd(tmp_dir)
clear tmp_dir
