function [g,filtertype] = gsp_jtv_design_damped_wave(G,alpha,beta,param)
%GSP_JTV_DESIGN_DAMPED_WAVE Design a damped wave time-vertex filter
%   Usage: [g,filtertype] = gsp_jtv_design_damped_wave(G);
%          [g,filtertype] = gsp_jtv_design_damped_wave(G,alpha);
%          [g,filtertype] = gsp_jtv_design_damped_wave(G,alpha,param);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       alpha   : Velocity parameters (default 1) 
%       beta    : Damping factors (default 0.1)
%       param   : Structure of optional parameters
%   Output parameters:
%       g          : Cell array of time-vertex filter
%       filtertype : Filter domain ts
%
%   gsp_jtv_design_damped_wave designs a damped wave filters according to the following eq. g(x,t) = exp(-beta abs(t))*cos(alpha t acos(1-x))
%   If alpha is a vector, a time-vertex filterbank will be designed. beta is a scalar.
%   Stability condition: alpha < 2/*fs*
%
%   Additional parameters
%   ---------------------
%    param.verbose*: verbosity level. 0 no log - 1 display warnings. (default 1)
%     
%   See also: gsp_jtv_design_wave, gsp_design_heat, gsp_jtv_design_dgw
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_jtv_design_damped_wave.php

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
% Date : July 2016

if nargin < 4
    param = struct;
end


if nargin<3
    if nargin<2
        alpha = 1;
        beta = 0.1;
    else
        beta = 0.1;
    end
end


if ~isfield(param,'verbose'), param.verbose = 1; end

K = gsp_jtv_design_wave(G,alpha,param);

H = gsp_design_heat(G.jtv.fs,beta);

g = gsp_jtv_design_dgw(G,K,@(x)1,H);

filtertype = 'ts';

end


