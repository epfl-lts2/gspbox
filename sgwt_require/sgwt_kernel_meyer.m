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
%
%   Url: http://lts2research.epfl.ch/gsp/doc/sgwt_require/sgwt_kernel_meyer.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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



