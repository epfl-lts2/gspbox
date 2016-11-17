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
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_norm_tik.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

