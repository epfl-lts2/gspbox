% sgwt_kernel_simple_tf : evaluates "simple" tight-frame kernel
%
% this is similar to meyer kernel, but simpler
%
% function is essentially sin^2(x) in ascending part,
% essentially cos^2 in descending part.
%
% function r= sgwt_kernel_simple_tf(x,kerneltype)
%
% Inputs
% x : array of independent variable values
% kerneltype : string, either 'sf' or 'wavelet' 
%
% Ouputs
% r : array of function values, same size as x.
%
% simple tf wavelet kernel : supported on [1/4,1]
% simple tf scaling function kernel : supported on [0,1/2]
%

function r= sgwt_kernel_simple_tf(x,kerneltype)
% h : [0,1]->[0,1] must satisfy h(0)=0, h(1)=1 .
h=@(x) sin(pi*x/2).^2;

r1ind=find(x>=0 & x<0.25);
r2ind=find(x>=.25 & x<0.5);
r3ind=find(x>=.5 & x<1);

r=zeros(size(x));

switch kerneltype
  case 'sf'
    r(r1ind)=1;
    r(r2ind)=sqrt(1-h(4*x(r2ind)-1).^2);
  case 'wavelet'
    r(r2ind)=h(4*(x(r2ind)-1/4));
    r(r3ind)=sqrt(1-h(2*x(r3ind)-1).^2);
end
