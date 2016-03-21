function [ gd ] = gsp_design_can_dual( g ,tol)
%GSP_DESIGN_CAN_DUAL This function return the canonical dual filters of g
%   Usage:  gd = gsp_design_can_dual( g );
%
%   Inputs parameters:
%       g       : cell array of filters
%       tol     : tolerance for the pseudo-inverse
%
%   Ouputs parameters:
%       g       : cell array of filters
%
%   This function returns the canonical dual filterbank g. Note that it
%   might not be the be the optimal solution in term of computation.
%
%   Example:
%
%             N = 100;
%             G = gsp_sensor(N);
%             G = gsp_compute_fourier_basis(G);
%             g = gsp_design_abspline(G,8);
%             gd = gsp_design_can_dual(g);
%             paramplot.show_sum = 0;
%             figure(1)
%             gsp_plot_filter(G,g,paramplot);
%             title('Original filters')
%             figure(2)
%             gsp_plot_filter(G,gd,paramplot);
%             title('Canonical dual filters');
% 
%             x = rand(N,1);
%             param.method = 'exact';
%             coeff = gsp_filter_analysis(G,g,x,param);
%             xs = gsp_filter_synthesis(G,gd,coeff,param);
%             norm(xs-x)
%
%   See also: gsp_evaluate_can_dual
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_can_dual.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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

if nargin<2
    tol = 1e-8;
end

Nf = length(g);
gd = cell(Nf,1);

for ii = 1:Nf
    gd{ii} = @(x) can_dual(g,ii,x,tol);
end
    
end


function ret = can_dual(g,n,x,tol)
    [N1, N2] = size(x);
    x = x(:);
    sol = gsp_evaluate_can_dual( g,x,tol );
    ret = sol(:,n);
    ret = reshape(ret,N1,N2);
end

