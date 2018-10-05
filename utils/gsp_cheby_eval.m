function r = gsp_cheby_eval(x,c,arange)
%GSP_CHEBY_EVAL Evaluate chebyshev polynomial
%   Usage:  r = gsp_cheby_eval(x,c,arange)
%
%   Input parameters:
%       x       : Points to evaluate the polynomial
%       c       : Chebyshef coefficients
%       arrange : arange (range to evaluate the polynomial)
%   Output parameters
%       r       : Result
%
%   In this function, arrange has to be [0, lmax ]!
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/utils/gsp_cheby_eval.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.5
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

% Author: David K Hammond, Nathanael Perraudin, Li Fan, David Shuman
% Testing: test_dual
% Date: 5 October 2018

% By setting the operator L to mean (pointwise) multiplication by x,
% and f to be vector of ones, p(L)f will contain p(x) at each
% point. This lets us use gsp_cheby_eval to evaluate the Chebyshev polynomials.


[N1,N2] = size(x);

L=spdiags(x(:),0,speye(numel(x)));
f=ones(size(x(:)));
N = length(f);

G.lmin = arange(1); 
G.lmax = arange(2);
G.L = L;
G.N = N;
r = gsp_cheby_op(G, c, f);

r = reshape(r, N1, N2);

end


