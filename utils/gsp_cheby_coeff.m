function c = gsp_cheby_coeff(G, filter, m, N,param)
%GSP_CHEBY_COEFF : Compute Chebyshev coefficients for a filterbank
%   Usage: c = gsp_cheby_coeff(G, filter, m, N);
%          c = gsp_cheby_coeff(G, filter, m);
%          c = gsp_cheby_coeff(G, filter);
%
%   Input parameters:
%       G       : graph structure or range of application
%       filter  : filter or cell array of filters
%       m       : maximum order Chebyshev coefficient to compute (default 30)
%       N       : grid order used to compute quadrature (default is m+1)
%       param   : optional parameter
%   Output parameters
%       c   : matrix of Chebyshev coefficients
% 
%   This function compute the Chebyshef coefficients for all the filter
%   contained in the cell array filter. The coefficient are returned in a
%   matrix. Every collumn correspond to a filter. The coefficients are
%   ordered such that c(j+1) is j'th Chebyshev coefficient
%
%   param contain only one field param.verbose to controle the verbosity.
%
%   Example:
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_meyer(G, Nf);  
%         c = gsp_cheby_coeff(G, g);
%
%   This function is inspired by the sgwt_toolbox
%
%   See also: gsp_cheby_op gsp_filter_analysis
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_cheby_coeff.php

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

% Author: David K Hammond, Nathanael Perraudin
% Testing: test_filter
% Date: 19 March 2014

if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 3
    m = 30;
end

if nargin < 4
   N = m+1; 
end

if iscell(filter)
   Nf = length(filter);
   c = zeros(m+1,Nf);
   for ii = 1: Nf
       c(:,ii) = gsp_cheby_coeff(G, filter{ii}, m, N,param);
   end
   return;
end

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_CHEBY_COEFF: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
        end
    end
  arange = [0, G.lmax];
else
  arange = G;
end
  
a1=(arange(2)-arange(1))/2;
a2=(arange(2)+arange(1))/2;
c = zeros(m+1,1);
for ii=1:m+1
    c(ii) = sum( filter( a1* cos( (pi*((1:N)-0.5))/N) + a2) .* ...
             cos( pi*(ii-1)*((1:N)-.5)/N) ) *2/N;
end

end
  

