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
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.normalize*   : Normalize the kernel
%
%

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

