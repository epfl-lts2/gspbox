function [ g,t ] = gsp_design_mexican_hat(G, Nf, param)
%GSP_DESIGN_MEXICAN_HAT Design the mexican hat filter bank
%   Usage: g =  gsp_design_mexican_hat(G, Nf, param);
%               gsp_design_mexican_hat(G ,Nf);
%               gsp_design_mexican_hat(G);
%
%   Input parameters:
%         G             : Graph or upper bound on the Laplacian spectrum
%         Nf            : Number of filters to cover the interval [0,lmax] (default 6)
%         param         : Structure of optional parameters
%   Output parameters:
%         g             : A cell array of filters
%
%   This function return a array of filters designed to be mexican hat
%   wavelet. The mexican hat wavelet is the second oder derivative of a
%   Gaussian. Since we express the filter in the Fourier domain, we find:
%
%           g_h(f) =   f^2 * exp(-f^2)
%
%      math: g_h(f) =  f^2 * exp(-f^2)
%
%   In our convention the eigenvalues of Laplacian are equivalent to the
%   square of vertex frequencies: f = lambda^2.
%
%   The low pass filter is given by
%
%           g_l(f) =   exp(-f^8)
%
%      math: g_l(f) =  exp(-f^8)
%
%   param is an optional structure containing the following fields
%
%    param.t*: vector of scale to be used (default: log scale)
%    param.lpfactor*: lmin*=*lmax*/*lpfactor will be used to determine
%     scales, then scaling function kernel will be created to fill the
%     lowpass gap. (default 20)
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
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
%         g = gsp_design_mexican_hat(G, Nf);   
%         gsp_plot_filter(G,g);  
%
%   This function is inspired by the sgwt_toolbox. 
%       
%   See also:
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_mexican_hat.php

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

% Author: Nathanael Perraudin, David K. Hammond
% Date: 18 March 2014

    

if nargin < 3
    param = struct;
end

if nargin < 2
    Nf = 6;
end

if ~isfield(param,'lpfactor'), param.lpfactor = 20; end
if ~isfield(param,'verbose'), param.verbose = 1; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_MEXICAN_HAT has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end

if ~isfield(param,'t')
   lmin=lmax / param.lpfactor;
   param.t = gsp_wlog_scales(lmin, lmax, Nf-1);
end


if param.verbose
    if length(param.t) ~= Nf - 1
       warning(['GSP_KERNEL_MEXICAN_HAT: You have specified ',...
           'more scales than Number of filter -1']);
    end
end

t = param.t;

% High pass filter
gb=@(x) x.*exp(-x);
% low pass filter
gl = @(x) exp(-x.^4);   

g = cell(Nf,1);

lminfac=.4*lmin;
g{1}=@(x) 1.2*exp(-1)*gl(x/lminfac);      

for j=1:Nf-1
    g{j+1}=@(x) gb(t(j)*x);
end



end


