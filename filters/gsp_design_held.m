function [ g ] = gsp_design_held(G, param)
%GSP_DESIGN_HELD Create a Held filterbank
%   Usage: g = gsp_design_held( G );
%          g = gsp_design_held( G, param );
%   
%   Inputs parameters:
%       G       : Graph structure or lmax
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbank
%
%   This function create a parseval filterbank of 2 filters. The low-pass
%   filter is defined by a function f_l(x): 
%
%                       /  1                                if x <= a 
%             f_l(x) = |   sin(2*pi*mu(val(r2ind)/(8*a)))   if a < x <= 2a
%                       \  0                                if x > 2a
%
%   with 
%
%             mu(x) = -1 + 24*x - 144*x^2 + 256*x^3
%   
%   The high pass filter is adaptated to obtain a tight frame.
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using:
%
%       G = gsp_estimate_lmax(G);
%
%   Example:
%
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_held(G);   
%         gsp_plot_filter(G,g);  
%         [A,B] = gsp_filterbank_bounds(G,g)
%
%   param is an optional structure containing the following fields
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%    param.a*: see equations above for this parameter. Note that the
%     spectrum is scaled between 0 and 2 (default 2/3).
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_held.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

% Author: Nathanael Perraudin, David Shuman
% Date  : 21 June 2014
% Testing: test_filter


if nargin < 2
    param = struct;
end


if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'a'), param.a = 2/3; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_HELD has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end




a = param.a;

g = cell(2,1);
g{1} = @(x) held(x*(2/lmax),a);
g{2} = @(x) real(sqrt(1-(held(x*(2/lmax),a)).^2));

end


function y = held(val,a)

y = zeros(size(val));

l1 = a;
l2 = 2*a;
mu = @(x) -1+24*x-144*x.^2+256*x.^3; 

r1ind = val >= 0     &    val < l1;
r2ind = val >= l1    &    val < l2;
r3ind = val >= l2;


y(r1ind) = 1;
y(r2ind) = sin(2*pi*mu(val(r2ind)/(8*a)));
y(r3ind) = 0;


end

