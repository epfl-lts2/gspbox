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
%   .. g(x) = s( (1-x) / bmax /lmax )
%
%   .. math:: g(x) = s\left(\frac{1-x}{\lambda_{\text{max}} b_{\text{max}} }\right)
%   
%   where  $s(x)$ is the step function
%
%   ..          /   0                                      if x < -1
%   ..  s(x) = | exp(-a/x) / ( exp(-a/x) + exp(-a/(1-x)) ) if x in [-1, 1]
%   ..          \   1                                      if x > 1
%
%   .. math:: s(x)=\begin{cases} 0 & \mbox{if }x<-1 \\ \frac{e^{-\frac{a}{x}}}{e^{-\frac{a}{x}}+e^{-\frac{a}{1-x}}} & \mbox{if }x\in[-1,1]\\ 1 & \mbox{if }x>1 \end{cases}
%
%   It uses a clever exponential construction to obtain an infinitely
%   differentiable function that is band limited!
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_expwin(G);   
%         gsp_plot_filter(G,g);  
%


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

g = @(x) gsp_smooth_downstep(x/bmax/lmax, a , 1);
    

    
end


