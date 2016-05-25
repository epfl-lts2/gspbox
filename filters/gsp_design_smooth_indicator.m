function [g]=gsp_design_smooth_indicator(G, a1, a2)
%GSP_DESIGN_SMOOTH_INDICATOR create a smooth indicator function
%   Usage: g = gsp_design_smooth_indicator(G, a1);
%          g = gsp_design_smooth_indicator(G, a1, a2);
%
%   Input parameters:
%       G       : Graph structure
%       a1      : Start of band cuttoff
%       a2      : End of band cuttoff (default 2*a1)
%
%   Output parameters
%       g       : filter
%
%   This function design the following filter:
%
%               /   1 if x < a1
%       g(x) = |    TBD
%               \   0 if x > a2
%
%   It use a clever exponential construction to obtain a infinitely
%   differentiable function that is band limited!
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
%         g = gsp_design_smooth_indicator(G,0.3,0.4);   
%         gsp_plot_filter(G,g);  
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_smooth_indicator.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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
% Date  : 18 December 2014


if nargin<2
    a1 = 0.25;
end

if nargin<3
    a2 = 2*a1;
end



if numel(a1)>1
    g = cell(numel(a1),1);
    for ii = 1:numel(a1);
        g{ii} = gsp_design_expwin(G,a1(ii),a2(ii));
    end
    return
end

if isstruct(G)
    if ~isfield(G,'lmax')
%         if param.verbose
            fprintf('GSP_DESIGN_HEAT has to compute lmax \n')
%         end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end

if a1>=a2
    error('GSP_DESIGN_SMOOTH_INDICATOR: a2 has to be bigger than a1!')
end

g = @(x) ffin(x/lmax,a1,a2);
    

    
end


function y=fx(x,a)  
    y = exp(-a./x);
    y(x<0)=0;
end

function y =gx(x,a)
    y=fx(x,a);
    y = y./(y+fx(1-x,a));
end

function y = ffin(x,a1,a2)
        y = gx(1-(x-a1)/(a2-a1),1);
        y(x<a1) = 1;
end
