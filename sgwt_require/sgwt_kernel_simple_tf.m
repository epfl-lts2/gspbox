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
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/sgwt_require/sgwt_kernel_simple_tf.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

function r= sgwt_kernel_simple_tf(x,kerneltype)
% h : [0,1]->[0,1] must satisfy h(0)=0, h(1)=1 .
h=@(x) sin(pi*x/2).^2;

%r1ind=find(x>=0 & x<0.25);
r1ind=find(x<0.25);

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

