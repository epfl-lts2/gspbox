function [ g,t ] = gsp_design_simple_tf(G, Nf, param)
%GSP_DESIGN_SIMPLE_TF Design a simple tight frame filterbank
%   Usage: g =  gsp_design_simple_tf(G, Nf, param);
%               gsp_design_simple_tf(G ,Nf);
%               gsp_design_simple_tf(G);
%
%   Input parameters:
%         G             : Graph or upper bound on the Laplacian spectrum
%         Nf            : Number of filters to cover the interval [0,lmax] (default 6)
%         param         : Structure of optional parameters
%   Output parameters:
%         g             : A cell array of filters
%
%   This function returns a array of filters designed to be simple tight
%   frame wavelet filterbank.
%
%   *param* is an optional structure containing the following fields
%
%   * *param.t*: vector of scale to be used (default: log scale)
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
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
%         g = gsp_design_simple_tf(G, Nf);   
%         gsp_plot_filter(G,g);  
%
%   This function is inspired by the sgwt_toolbox. 
%       
%   See also:

% Author: Nathanael Perraudin, David K. Hammond
% Date: 18 March 2014

    

if nargin < 3
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_SIMPLE_TF has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end

if ~isfield(param,'t')
   param.t = (1/(2*lmax)) * 2.^(Nf-2:-1:0);
end


if param.verbose
    if length(param.t) ~= Nf - 1
       warning(['GSP_KERNEL_SIMPLE_TF: You have specified ',...
           'more scales than Number of filter -1']);
    end
end

t = param.t;
g = cell(Nf,1);


g{1}= @(x) sgwt_kernel_simple_tf(t(1)*x,'sf');
for j=1:Nf-1
    g{j+1}= @(x) sgwt_kernel_simple_tf(t(j)*x,'wavelet');
end



end

