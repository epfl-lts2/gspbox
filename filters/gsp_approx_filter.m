function gt = gsp_approx_filter(G, g, m, N, param)
%GSP_APPROX_FILTER : Create the approximation filter for a filterbank
%   Usage: c = gsp_approx_filter(G, filter, m, N);
%          c = gsp_approx_filter(G, filter, m);
%          c = gsp_approx_filter(G, filter);
%          c = gsp_approx_filter(G, filter, m, N,param);
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
%   This function create the approximate filters of g with a truncated
%   Chebyshev polynomial.
%
%   param contain only one field param.verbose to controle the verbosity.
%
%   Example:
%
%             N = 100;
%             order = 15;
%             G = gsp_sensor(N);
%             G = gsp_estimate_lmax(G);
%             g = gsp_design_abspline(G,8);
%             ga = gsp_approx_filter(G,g,order);
%             paramplot.show_sum = 0;
%             figure(1)
%             gsp_plot_filter(G,g,paramplot);
%             title('Original filters')
%             figure(2)
%             gsp_plot_filter(G,ga,paramplot);
%             title('Approximate filters');
%             x = rand(N,1);
%             param.order = order;
%             c1 = gsp_filter_analysis(G,g,x,param);
%             c2 = gsp_filter_analysis(G,ga,x,param);
%             norm(c1-c2)/norm(c1)
%
%   See also: gsp_cheby_eval gsp_filter_analysis
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_approx_filter.php

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
% Date  : 30 December 2014
% Testing: test_dual

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

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_APPROX_FILTER: The variable lmax is not ',...
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
  

Nf = length(g);

c = gsp_cheby_coeff(G, g, m, N,param);

gt = cell(Nf,1);

for ii = 1:Nf
    gt{ii} = @(x) gsp_cheby_eval(x,c(:,ii),arange);
end

if Nf == 1;
    gt = gt{1};
end

end

