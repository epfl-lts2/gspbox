function [ g,mu ] = gsp_design_itersine( G, Nf, param )
%GSP_DESIGN_ITERSINE Create a itersine filterbanks
%   Usage: g = gsp_design_itersine( G, Nf );
%          g = gsp_design_itersine( G, Nf, param );
%   
%   Inputs parameters:
%       G       : Graph structure or lmax
%       Nf      : Number of filter 
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbanks
%       mu      : centers of the filters
%
%   This function creates an itersine half overlap filterbank of *Nf* filters
%   Going from $0$ to $\lambda_{max}$
%
%   This filterbank is tight for an overlap of 2 and other particular
%   values. The function normalizes the window such that the framebound is
%   1.
%
%   The itersine window between -0.5 and 0.5 is defined as
%
%   ..  g(t) = sin( 0.5 pi (cos pi t)^2 )
%
%   .. math:: g(t) = \sin \left( 0.5 \pi \cos (\pi t)^2 \right) 
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%         Nf = 20;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_itersine(G, Nf);   
%         gsp_plot_filter(G,g);  
%         [A,B] = gsp_filterbank_bounds(G,g)
%
%   *param* is an optional structure containing the following fields
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.overlap*: Overlap : default 2
%

% Author: Nathanael Perraudin
% Date  : 18 June 2014
% Testing: test_filter

if nargin < 3
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'overlap'), param.overlap = 2; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_ITERSINE has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end


k = @(x) sin(0.5 * pi * (cos(pi*x)).^2) .* (x>=-0.5 & x<= 0.5);

g = cell(Nf,1);

scale = lmax/(Nf-param.overlap+1)*(param.overlap);

mu = zeros(Nf,1);

for ii=1:Nf
    g{ii} = @(x) k(x/scale-(ii-param.overlap/2)/param.overlap)...
                    ./sqrt(param.overlap)*sqrt(2);
    mu(ii) = (ii-param.overlap/2)/param.overlap * scale;
end


end

