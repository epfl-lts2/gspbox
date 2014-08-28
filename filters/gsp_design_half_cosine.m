function [ filters ] = gsp_design_half_cosine( G, Nf,param )
%GSP_DESIGN_HALF_COSINE Design uniform half cosine filterbank
%   Usage: filters = gsp_design_half_cosine( Nf, UBT );
%
%   Inputs parameters:
%       G       : Graph or maximum value
%       Nf      : Number of filters
%
%   Outputs parameters:
%       filters : Cell array of filters
%
%   This function generate a uniform half cosine filterbank. The main
%   window
%
%       0.5 * (1+cos(2*pi*(x/a-1/2)))  for  0 <= x <= a
%
%   is translated uniformly to create the filterbank.
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using:
%
%       G = gsp_estimate_lmax(G);
%
%   Example:
% 
%         figure(100);
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_half_cosine(G, Nf);   
%         gsp_plot_filter(G,g);
%
%   param is an optional structure containing the following fields
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_half_cosine.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

% Author: David Shuman, Nathanael Perraudin
% Date  : 15 June 2014
% Testing: test_filter


if nargin < 3
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_HALF_COSINE has to compute lmax \n')
        end
            G = gsp_estimate_lmax(G);
    end
    lmax = G.lmax;
else
    lmax = G;
end

% Design main window
dilation_factor = lmax*(3/(Nf-2));
main_window = @(x) (.5+.5*cos(2*pi*(x/dilation_factor-1/2)))...
                   .*(x>=0).*(x<=dilation_factor);

% Design filters (uniform translates of the main window, cut-off at the
% spectrum boundaries) 

filters = cell(Nf,1);
for i=1:Nf
   filters{i} = @(x) main_window(x-dilation_factor/3*(i-3)); 
end

end


