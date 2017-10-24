function [g,filtertype] = gsp_jtv_design_diffusion(G,tau,param)
%GSP_JTV_DESIGN_DIFFUSION Design a diffusion time-vertex filter
%   Usage: [g,ft] = gsp_jtv_design_diffusion(G);
%          [g,ft] = gsp_jtv_design_diffusion(G,tau);
%          [g,ft] = gsp_jtv_design_diffusion(G,tau,param);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       tau     : Diffusivity parameter (default 1)
%       param   : Struct of optional paramenter
%
%   Output parameters
%       g          : Cell array of time-vertex filters
%       filtertype : Filter domain: ts
%
%   Design a diffusion time-vertex filter
%
%   Additional parameters
%   ---------------------
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%    param.normalize   : Normalize the kernel
%
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_design_diffusion.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

% Author :  Francesco Grassi
% Date : July 2016


if nargin < 3
    param = struct;
end

if nargin<2
    tau=1;
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'normalize'), param.normalize = 0; end

if numel(tau)>1
    Nf = numel(tau);
    g = cell(Nf,1);
    for ii = 1:Nf
        [g{ii},filtertype] = gsp_jtv_design_diffusion(G,tau(ii),param); 
    end
    return
end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            warning('GSP_JTV_DESIGN_DIFFUSION has to compute lmax')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
   if ~gsp_check_jtv(G)
        error('GSP_JTV_DESIGN_DIFFUSION need the time dimension. Use GSP_JTV_GRAPH');
    end
   T = G.jtv.T;
   fs = G.jtv.fs;
   
else
   lmax = G;
end

if param.normalize
    
    filter_norm = -psi(1) + log(2*tau*T/fs) + expint(2*tau*T/fs);
    
    g = @(x,t) (t<(T/fs)).*(t>=0).*exp(-tau*abs(t).*x/(lmax*fs))/filter_norm;
    
else
    
    g = @(x,t) (t<(T/fs)).*(t>=0).*exp(-tau*abs(t).*x/(lmax*fs));
    
end

filtertype = 'ts';

end


