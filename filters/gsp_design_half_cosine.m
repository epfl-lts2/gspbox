function [ filters ] = gsp_design_half_cosine( G, Nf,param )
%GSP_DESIGN_HALF_COSINE Design uniform half cosine filterbank
%   Usage: filters = gsp_design_half_cosine( Nf, UBT );
%
%   Inputs parameters:
%       G       : Graph or maximum value
%       Nf      : Number of filters
%
%   Outputs parameters:
%       filters : Cell array of filters
%
%   This function generates a uniform half cosine filterbank. The main
%   window
%
%   ..  0.5 * (1+cos(2*pi*(x/a-1/2)))  for  0 <= x <= a
%
%   .. math:: \frac{1}{2} \left(1 + \cos\left(2\pi\left(\frac{x}{a}-\frac{1}{2}\right)\right)\right)  \text{for } 0 \leq x \leq a
%
%   is translated uniformaly to create the filterbank.
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
%         g = gsp_design_half_cosine(G, Nf);   
%         gsp_plot_filter(G,g); 
%
%   *param* is an optional structure containing the following fields
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%

% Author: David Shuman, Nathanael Perraudin
% Date  : 15 June 2014
% Testing: test_filter


if nargin < 3
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_HALF_COSINE has to compute lmax \n')
        end
            G = gsp_estimate_lmax(G);
    end
    lmax = G.lmax;
else
    lmax = G;
end

% Design main window
dilation_factor = lmax*(3/(Nf-2));
main_window = @(x) (.5+.5*cos(2*pi*(x/dilation_factor-1/2)))...
                   .*(x>=0).*(x<=dilation_factor);

% Design filters (uniform translates of the main window, cut-off at the
% spectrum boundaries) 

filters = cell(Nf,1);
for i=1:Nf
   filters{i} = @(x) main_window(x-dilation_factor/3*(i-3)); 
end

end

