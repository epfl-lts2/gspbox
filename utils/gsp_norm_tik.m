function y = gsp_norm_tik(G,x)
%GSP_NORM_TIK Squared L2 norm of the gradient on graph
%   Usage:  y = gsp_norm_tv(G,x);
%
%   Input parameters:
%         G     : Graph structure (or symetric positive matrix)
%         x     : Signal on graph
%   Output parameters:
%         y     : Norm
%
%   Compute the squared L2 norm of the gradient on graph. If x is a matrix
%   a vector of norm is returned.
%
%   This function can also be used for general symetric positive matrices
%
%   See also: gsp_prox_tik gsp_norm_tv

% Author: Nathanael Perraudin
% Date:   25 March 2014

if isa(x,'single')
   x = double(x); 
end

if ~isnumeric(G)
    L = G.L;
else
    L = G;
end



[N,M ] = size(L);
if size(x,1) ~= M
    error('The dimension of x is not compatible with the dimension of L');
end
NL = M/N;
y = 0;
for ii = 1:NL;
    ind = (1:N)+(ii-1)*N;
    y = y + sum(x(ind,:) .* (L(:,ind)* x(ind,:)) );
end

% Previous implementation
% y = sum(x .* (L* x) );

end
