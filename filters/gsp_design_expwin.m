function [g]=gsp_design_expwin(G, bmax, a)
%GSP_DESIGN_EXPWIN create an expwin window of lenth N with parameter a
%   Usage: g = gsp_design_expwin(G);
%          g = gsp_design_expwin(G, bmax);
%          g = gsp_design_expwin(G, bmax, a);
%
%   Input parameters:
%       G       : Graph structure
%       bmax    : Maximum relative band (default 0.2)
%       a       : Slope parameter (default 1).
%
%   Output parameters
%       g       : filter
%
%   This function design the following filter:
%
%       g(x) = 
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
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_expwin(G);   
%         gsp_plot_filter(G,g);  
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_expwin.php

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
% Date  : 18 December 2014


if nargin<3
    a = 1;
end


if nargin<2
    bmax = 0.2;
end


if numel(bmax)>1
    g = cell(numel(bmax),1);
    for ii = 1:numel(bmax);
        g{ii} = gsp_design_expwin(G,bmax(ii),a);
    end
    return
end

if isstruct(G)
    if ~isfield(G,'lmax')
%         if param.verbose
            fprintf('GSP_DESIGN_EXPWIN has to compute lmax \n')
%         end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end


g = @(x) ffin(x/bmax/lmax,a);
    

    
end


function y=fx(x,a)
    y = exp(-a./x);
    y(x<0)=0;
end

function y =gx(x,a)
    y=fx(x,a);
    y = y./(y+fx(1-x,a));
end

function y = ffin(x,a)
        y = gx(1-x,a);
end

