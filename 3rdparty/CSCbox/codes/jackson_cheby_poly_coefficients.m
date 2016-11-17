% function  [CH,JCH] = jackson_cheby_poly_coefficients(a,b,lambda_range,m)
%
% That's to compute the m+1 coefficients of the polynomial approximation of 
% an ideal band-pass between a and b, in between a range of values defined by 
% lambda_range=[lambda_min,lambda_max]; 
% ----
% Output:
% - CH are the coefficients of the Chebychev polynomials
% - JCH are the coefficients of the jackson-chebychev polynomials
% ----
% For details of the following calculations, see 
% L. O. Jay et al., "Electronic structure calculations for plane-wave codes without diagonalization",
% published in "Computer Physics Communications ", vol. 118, no. 1, pp. 21-30, 1999.
%
% Copyright (C) 2016 Nicolas Tremblay, Gilles Puy.
% This file is part of the CSCbox (Compressive Spectral Clustering toolbox)
%
% The CSCbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The CSCbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%     N. Tremblay, G. Puy, R. Gribonval and P. Vandergheynst.
%     Compressive Spectral Clustering.
%     ArXiv e-prints:1602.02018 Feb. 2016.

function  [CH,JCH] = jackson_cheby_poly_coefficients(a,b,lambda_range,m)

% scaling and translation coefficients compared to the classical interval
% of Chebychev polynomials [-1,1] :
a1 = (lambda_range(2)-lambda_range(1))/2;
a2 = (lambda_range(1)+lambda_range(2))/2;

% scale the boundaries of the band pass according to lrange:
a=(a-a2)/a1;
b=(b-a2)/a1;

% compute Cheby coef:
CH(1)=(1/pi)*(acos(a)-acos(b));
for j=2:m+1
    CH(j)=(2/(pi*(j-1)))*(sin((j-1)*acos(a))-sin((j-1)*acos(b)));
end

% compute Jackson coef:
alpha=pi/(m+2);
for j=1:m+1
    gamma_JACK(j)=(1/sin(alpha))*((1-(j-1)/(m+2))*sin(alpha)*cos((j-1)*alpha)+(1/(m+2))*cos(alpha)*sin((j-1)*alpha));
end

% compute Jackson-Cheby coef:
JCH=CH.*gamma_JACK;

% to be in adequation with gsp_cheby_op.m :
JCH(1)=JCH(1)*2;
CH(1)=CH(1)*2;

JCH=JCH';
CH=CH';
