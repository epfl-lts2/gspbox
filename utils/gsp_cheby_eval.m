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
%   See also: gsp_cheby_op
%

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


