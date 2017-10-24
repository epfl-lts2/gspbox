% sgwt_kernel_abspline5 : Monic polynomial / quintic spline / power law decay kernel
%
% function r = sgwt_kernel_abspline5(x,alpha,beta,t1,t2)
%
% Defines function g(x) with g(x) = c1*x^alpha for 0<x<x1
% g(x) = c3/x^beta for x>t2
% quintic spline for t1<x<t2,
% Satisfying g(t1)=g(t2)=1
% g'(t1)=g'(t2)
% g''(t1)=g''(t2)
%
% Inputs :
% x : array of independent variable values
% alpha : exponent for region near origin
% beta : exponent decay
% t1, t2 : determine transition region
%
% Outputs :
% r - result (same size as x)

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function r = sgwt_kernel_abspline5(x,alpha,beta,t1,t2)
  r=zeros(size(x));
  % compute spline coefficients
  % M a = v
  M=[[1 t1 t1^2 t1^3 t1^4 t1^5];...
     [1 t2 t2^2 t2^3 t2^4 t2^5];...
     [0 1 2*t1 3*t1^2 4*t1^3 5*t1^4];...
     [0 1 2*t2 3*t2^2 4*t2^3 5*t2^4];
     [0 0 2 6*t1 12*t1^2 20*t1^3];...
     [0 0 2 6*t2 12*t2^2 20*t2^3]...
     ];
  %v=[t1^alpha ; t2^(-beta) ; alpha*t1^(alpha-1) ; -beta*t2^(-beta-1)];
  v=[1 ; 1 ; ...
     t1^(-alpha)*alpha*t1^(alpha-1) ; -beta*t2^(-beta-1)*t2^beta; ...
     t1^(-alpha)*alpha*(alpha-1)*t1^(alpha-2);-beta*(-beta-1)*t2^(-beta-2)*t2^beta
     ];
  a=M\v;
  
  r1=find(x<t1);
  %r1=find(x>=0 & x<t1);

  r2=find(x>=t1 & x<t2);
  r3=find(x>=t2);
  r(r1)=x(r1).^alpha*t1^(-alpha);
  r(r3)=x(r3).^(-beta)*t2^(beta);
  
  x2=x(r2);
  r(r2)=a(1)+a(2)*x2+a(3)*x2.^2+a(4)*x2.^3+a(5)*x2.^4+a(6)*x2.^5;
%  tmp=polyval(flipud(a),x2);
%  keyboard
