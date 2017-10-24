% sgwt_kernel_meyer : evaluates meyer wavelet kernel and scaling function
% function r=sgwt_kernel_meyer(x,kerneltype)
%
% Inputs
% x : array of independent variable values
% kerneltype : string, either 'sf' or 'wavelet' 
%
% Ouputs
% r : array of function values, same size as x.
%
% meyer wavelet kernel : supported on [2/3,8/3]
% meyer scaling function kernel : supported on [0,4/3]
%
% Use of this kernel for SGWT proposed by Nora Leonardi and Dimitri Van De Ville,
% "Wavelet Frames on Graphs Defined by fMRI Functional Connectivity"
% International Symposium on Biomedical Imaging, 2011
function r=sgwt_kernel_meyer(x,kerneltype)
l1=2/3;
l2=4/3;%2*l1;
l3=8/3;%4*l1;
v=@(x) x.^4.*(35-84*x+70*x.^2-20*x.^3) ; 

r1ind=find(x>=0 & x<l1);
r2ind=find(x>=l1 & x<l2);
r3ind=find(x>=l2 & x<l3);
% as we initialize r with zero, computed function will implicitly be zero for
% all x not in one of the three regions defined above
r=zeros(size(x));
switch kerneltype
  case 'sf'
    r(r1ind)=1;
    r(r2ind)=cos((pi/2)*v(abs(x(r2ind))/l1-1));
  case 'wavelet'
    r(r2ind)=sin((pi/2)*v(abs(x(r2ind))/l1-1));
    r(r3ind)=cos((pi/2)*v(abs(x(r3ind))/l2-1));
  otherwise
    error(sprintf('unknown kernel type %s',kerneltype));
end


