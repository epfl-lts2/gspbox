function [ bool ] = gsp_test_duality_coefficient( gcoeff,hcoeff,tol )
%GSP_TEST_DUALITY_COEFFICIENT Test if the coefficient are from dual filters
%   Usage:  bool = gsp_test_duality_coefficient( gcoeff,hcoeff );
%           bool = gsp_test_duality_coefficient( gcoeff,hcoeff,tol );
%
%   Input parameters:
%       gcoeff  : coefficient of the filter 1 (matrix $N$ x $M$ )
%       hcoeff  : coefficinet of the filter 2 (matrix $N$ x $M$ )
%       tol     : tolerance for the test (default 1e-5)
%
%   Ouput paramters:
%       bool    : boolean 
%
%   This function test if two discrete filterbanks are dual. Each filter is
%   a column in the matrix *gcoeff* or *hcoeff*. $M$ is the number of
%   filters and $N$ the number of coefficients (size of the graph signal).
%

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

