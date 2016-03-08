function [ g ] = gsp_design_regular(G, param)
%GSP_DESIGN_REGULAR Create a Simoncelli filterbank
%   Usage: g = gsp_design_regular( G );
%          g = gsp_design_regular( G, param );
%   
%   Inputs parameters:
%       G       : Graph structure or lmax
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbank
%
%   This function create a parseval filterbank of 2 filters. The low-pass
%   filter is defined by a function f_l(x) between 0 and 2. For
%   d = 0.
%
%         f_l(x) = sin(pi/4*x)
%
%   For d = 1 
%
%         f_l(x) = sin( pi/4 * (1+sin(pi/2*(x-1))) )
%
%   For d = 2 
%
%         f_l(x) = sin( pi/4 * ( 1 + sin( pi/2 * sin(pi/2*(x-1) ) ) )
%
%   And so on for the other degrees d.
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
%         g = gsp_design_regular(G);   
%         gsp_plot_filter(G,g);  
%         [A,B] = gsp_filterbank_bounds(G,g)
%
%   param is an optional structure containing the following fields
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%    param.d*: Degree. See equation for mor informations. (default 3)
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_regular.php

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

% Author: Nathanael Perraudin, David Shuman
% Date  : 21 June 2014
% Testing: test_filter


if nargin < 2
    param = struct;
end


if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'d'), param.d = 3; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_REGULAR has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end




d = param.d;

g = cell(2,1);
g{1} = @(x) regular(x*(2/lmax),d);
g{2} = @(x) real(sqrt(1-(regular(x*(2/lmax),d)).^2));

end


function y = regular(val,d)


if d==0
    y = sin(pi/4*val);
else
    output = sin(pi*(val-1)/2);
    for k=2:d
        output = sin(pi*output/2);
    end
    y = sin(pi/4*(1+output));
end


end

