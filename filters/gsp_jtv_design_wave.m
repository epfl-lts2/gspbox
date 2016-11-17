function [g,filtertype] = gsp_jtv_design_wave(G,alpha,param)
%GSP_JTV_DESIGN_WAVE Design a jtv wave filterbank
%   Usage: gsp_jtv_design_wave(G);
%          gsp_jtv_design_wave(G,alpha);
%          gsp_jtv_design_wave(G,alpha,param);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       alpha   : Velocity parameter (default 1)
%   Output parameters:
%       g          : Cell array of time-vertex filters
%       filtertype : Filter domain ts
%
%   Stability condition: alpha < 2/*fs*
%
%   Additional parameters
%   ---------------------
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1)
%    param.normalize   : Normalize the kernel
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_jtv_design_wave.php

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

% Author :  Francesco Grassi
% Date   : September 2016


if nargin < 3
    param = struct;
end

if nargin<2
    alpha=1;
end

if ~isfield(param,'normalize'), param.normalize = 0; end
if ~isfield(param,'verbose'), param.verbose = 0; end

if isstruct(G)
    if ~isfield(G.jtv,'T')
        error('GSP_JTV_DESIGN_WAVE needs the time dimension. Use GSP_JTV_GRAPH');
    end
    
    T = G.jtv.T;
    fs = G.jtv.fs;
    
    if ~isfield(G,'lmax')
        if param.verbose
            warning('GSP_JTV_DESIGN_WAVE has to compute lmax')
        end
        G = gsp_estimate_lmax(G);
    end
    
    lmax = G.lmax;
    
    CFL = 2/fs;
    
else
    error('Graph structure not valid');
end

if any(alpha>CFL)
    error('The filter is unstable. Try reducing alpha.');
end

if numel(alpha)>1
    Nf = numel(alpha);
    g = cell(Nf,1);
    for ii = 1:Nf
        [g{ii},filtertype] = gsp_jtv_design_wave(G,alpha(ii),param);
    end
    return
end


if param.normalize

    
    if alpha == 0
        fie_norm = T*G.lmax;
    else
        s = linspace( (2*T-1),(2*T+1),100)*acos(1-alpha^2/2);
        fie_norm = 2*G.lmax*(2*alpha^2*T + log((2*T+1)./(2*T-1)) - trapz(s,cos(s-1)./s))/(8*alpha^2);
    end
    
    g = @(x,t) (t<T).*(t>=0).*cos(t.*acos(1-alpha^2*x/(2*lmax)))/fie_norm;
    
else
    
    g = @(x,t) (t<T).*(t>=0).*cos(t.*acos(1-alpha^2*x/(2*lmax)));
    
end
filtertype = 'ts';

end


