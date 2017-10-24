function [ g,t ] = gsp_design_mexican_hat(G, Nf, param)
%GSP_DESIGN_MEXICAN_HAT Design the mexican hat filterbank
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
%   ..      g_h(f) =   f^2 * exp(-f^2)
%
%   .. math: g_h(f) =  f^2 * exp(-f^2)
%
%   In our convention the eigenvalues of Laplacian are equivalent to the
%   square of vertex frequencies: $f = \lambda^2$.
%
%   The low pass filter is given by
%
%   ..      g_l(f) =   exp(-f^8)
%
%   .. math: g_l(f) =  exp(-f^8)
%
%   *param* is an optional structure containing the following fields
%
%   * *param.t*: vector of scale to be used (default: log scale)
%   * *param.lpfactor*: *lmin*=*lmax*/*lpfactor* will be used to determine
%     scales, then scaling function kernel will be created to fill the
%     lowpass gap. (default 20)
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.normalize*: normalize the wavelet by the factor $\sqrt{t}$
%     (default 0.)
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%         Nf = 8;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_mexican_hat(G, Nf);   
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

if nargin < 2
    Nf = 6;
end

if ~isfield(param,'lpfactor'), param.lpfactor = 20; end
if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'normalize'), param.normalize = 0; end

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

lmin=lmax / param.lpfactor;

if ~isfield(param,'t')
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
    if param.normalize
        g{j+1} = @(x) sqrt(t(j))*gb(t(j)*x);
    else
        g{j+1} = @(x) gb(t(j)*x);
    end
end



end

