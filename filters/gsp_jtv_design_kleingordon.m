function [g,filtertype] = gsp_jtv_design_kleingordon(G,alpha,mu,param)
%GSP_JTV_DESIGN_KLEINGORDON Design a jtv Klein-Gordon filterbank
%   Usage: gsp_jtv_design_kleingordon(G);
%          gsp_jtv_design_kleingordon(G,alpha,mu);
%          gsp_jtv_design_kleingordon(G,alpha,mu,param);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       alpha   : Velocity parameter (default 1)
%   Output parameters:
%       g          : Cell array of time-vertex filters
%       filtertype : Filter domain ts
%
%
%   Additional parameters
%   ---------------------
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1)
%   * *param.normalize   : Normalize the kernel

% Author :  Francesco Grassi
% Date   : September 2016


if nargin < 4
    param = struct;
end

if nargin<3
    mu = .1;
end

if nargin<2
    alpha = 1;
end

if ~isfield(param,'normalize'), param.normalize = 0; end
if ~isfield(param,'verbose'), param.verbose = 0; end

if isstruct(G)
    if ~isfield(G.jtv,'T')
        error('GSP_JTV_DESIGN_KLEINGORDON needs the time dimension. Use GSP_JTV_GRAPH');
    end
    
    T = G.jtv.T;
    fs = G.jtv.fs;
    
    if ~isfield(G,'lmax')
        if param.verbose
            warning('GSP_JTV_DESIGN_KLEINGORDON has to compute lmax')
        end
        G = gsp_estimate_lmax(G);
    end
    
    lmax = G.lmax;
    
    CFL = 2/fs;
    
else
    error('Graph structure not valid');
end

% if any(alpha>CFL)
%     error('The filter is unstable. Try reducing alpha.');
% end

% if numel(alpha)>1
%     Nf = numel(alpha);
%     g = cell(Nf,1);
%     for ii = 1:Nf
%         [g{ii},filtertype] = gsp_design_wave(G,alpha(ii),param);
%     end
%     return
% end



g = @(x,t) (t<T).*(t>=0).*cos(t.*acos(1-alpha^2*(x+mu)/(2*lmax)));
    

filtertype = 'ts';

end

