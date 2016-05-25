function [ bool ] = gsp_test_duality_coefficient( gcoeff,hcoeff,tol )
%GSP_TEST_DUALITY_COEFFICIENT Test if the coefficient are from dual filters
%   Usage:  bool = gsp_test_duality_coefficient( gcoeff,hcoeff );
%           bool = gsp_test_duality_coefficient( gcoeff,hcoeff,tol );
%
%   Input parameters:
%       gcoeff  : coefficient of the filter 1 (matrix N x M )
%       hcoeff  : coefficinet of the filter 2 (matrix N x M )
%       tol     : tolerance for the test (default 1e-5)
%
%   Ouput paramters:
%       bool    : boolean 
%
%   This function test if two discrete filterbanks are dual. Each filter is
%   a culumn in the matrix gcoeff or hcoeff. M is the number of
%   filters and N the number of coefficients (size of the graph signal).
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_test_duality_coefficient.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

% Author: Nathanael Perraudin
% Date  : 13 july 2014
% testing: test_dual

if nargin<3
    tol = 1e-5;
end

v = sum(gcoeff.*hcoeff,2);
A = min(v);
B = max(v);
bool = ( (A-B) / ((A+B)/2) ) < tol;

end


