function [ y ] = gsp_frkhs( x,alpha,xi,k,method )
%GSP_FRKHS This function evalute the function f
%   Usage: y = gsp_frkhs( x,alpha,xi,k )
%
%   Input parameters
%       x       : column vector of point to evaluate the function
%       alpha   : coefficient of the function (vector column!)
%       xi      : measurements
%       k       : kernel
%       method  : method for the function rkhs_evalute (default mixed)
%   Output parameters
%       y       : column vector of solution
%
%   This function evaluates the function f in this form
%
%        y  = sum_l alpha_l k(x_l,x)
%
%   This is a dummy function that need improvements
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/gsp_frkhs.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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
% Date  : 10 june 2014

if nargin < 5
    method = 'mixed';
end

M = size(x,2);

y = zeros(M,size(alpha,2));
N = length(alpha);

if size(alpha,2)>1
    for ii = 1: size(alpha,2)
        y(:,ii) = gsp_frkhs( x,alpha(:,ii),xi,k);
    end
    return
end

if iscell(k)
    Nc = length(k);
    for ii = 1:Nc
            y = y + gsp_frkhs( x,alpha((1:(N/Nc))+(ii-1)*N/Nc),xi,k{ii} );
    end
    return
end


% 
% for ii = 1:N
%     y = y + alpha(ii) * k(repmat(xi(ii,:),M,1),x);
% end

Kx = gsp_rkhs_evaluate(k,xi,x,method);
y =   Kx'*alpha;

end


