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
%   Stability condition: *alpha* < 2/*fs*
%
%   Additional parameters
%   ---------------------
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings. (default 1)
%     
%   See also: gsp_jtv_design_wave, gsp_design_heat, gsp_jtv_design_dgw

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

