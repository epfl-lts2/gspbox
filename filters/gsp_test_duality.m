function [ bool ] = gsp_test_duality(G, g,h,tol )
%GSP_TEST_DUALITY Test if two filterbanks are dual
%   Usage: bool = gsp_test_duality(G, g,h )
%          bool = gsp_test_duality(G, g,h,tol )
%
%   Input parameters:
%       G       : Graph or arange (min and max value) 
%       g       : filter 1 (or filterbank)
%       h       : filter 2 (or filterbank)
%       tol     : tolerance for the test (default 1e-8)
%
%   Ouput paramters:
%       bool    : boolean 
%
%   This function test if two filterbanks are dual.
%
%   Example:
%
%             N = 100;
%             G = gsp_sensor(N);
%             G = gsp_estimate_lmax(G);
%             g = gsp_design_abspline(G,8);
%             gd = gsp_design_can_dual(g);
%             gsp_test_duality(G, g,gd )
%
%   See also: gsp_design_can_dual gsp_test_duality_coefficient
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_test_duality.php

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
% Date  : 30 December 2014
% Testing: test_dual

if nargin<4
    tol = 1e-8;
end

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_TEST_DUALITY: The variable lmax is not ',...
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

N = 100;

x = linspace(arange(1),arange(2),N)';
c1 = gsp_filter_evaluate(g,x);
c2 = gsp_filter_evaluate(h,x);

bool = gsp_test_duality_coefficient(c1,c2,tol);


end


