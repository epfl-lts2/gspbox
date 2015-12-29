function [ s ] = gsp_evaluate_can_dual( g,val,tol )
%GSP_EVALUATE_CAN_DUAL Evaluate the canonical dual filterbank
%   Usage: hcoeff = gsp_evaluate_can_dual( g,val )
%
%   Inputs parameters:
%       g       : cell array of filters
%       val     : column vectors of values
%       tol     : tolerance
%
%   Ouputs parameters:
%       s       : Matrix of value
%
%   This function compute the value of the canonical dual of a filterbank
%   g at the point specified in val. The function return a matrix. Each
%   column is the output of one dual filter.
%
%   See also: gsp_design_can_dual
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_evaluate_can_dual.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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
% Date  : 14 June 2014
% Testing: test_dual

if nargin<3
    tol = 1e-10;
end


% TODO: size should be improved
N = length(val);

% Compute coefficient of g
gcoeff = gsp_filter_evaluate(g,val)';

% Compute coefficient of h
% s = zeros(N,M);
% for ii = 1:N
%     s(ii,:) =  pinv(gcoeff(ii,:)'); 
% end

s = arrayfun(@(x) pinv(gcoeff(:, x),tol), 1 : N, 'UniformOutput', false);
s = cell2mat(s');
end


