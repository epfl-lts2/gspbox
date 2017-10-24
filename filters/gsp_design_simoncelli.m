function [ g ] = gsp_design_simoncelli(G, param)
%GSP_DESIGN_SIMONCELLI Create a Simoncelli filterbank
%   Usage: g = gsp_design_simoncelli( G );
%          g = gsp_design_simoncelli( G, param );
%   
%   Inputs parameters:
%       G       : Graph structure or lmax
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbank
%
%   This function creates a parseval filterbank of $2$ filters. The low-pass
%   filter is defined by a function $f_l(x)$: 
%
%   ..              /  1                                   if x <= a 
%   ..    f_l(x) = |   cos( pi/2 * log(x/a) / log(2) )     if a < x <= 2a
%   ..              \  0                                   if x > 2a
%
%   .. math:: f_{l}=\begin{cases} 1 & \mbox{if }x\leq a\\ \cos\left(\frac{\pi}{2}\frac{\log\left(\frac{x}{2}\right)}{\log(2)}\right) & \mbox{if }a<x\leq2a\\ 0 & \mbox{if }x>2a \end{cases}
%
%   The high pass filter is adaptated to obtain a tight frame.
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%
%   Example:::
%
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_simoncelli(G);   
%         gsp_plot_filter(G,g);  
%         [A,B] = gsp_filterbank_bounds(G,g)
%
%   *param* is an optional structure containing the following fields
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.a*: see equations above for this parameter. Note that the
%     spectrum is scaled between 0 and 2 (default 2/3).
%

% Author: Nathanael Perraudin, David Shuman
% Date  : 21 June 2014
% Testing: test_filter


if nargin < 2
    param = struct;
end


if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'a'), param.a = 2/3; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_SIMONCELLI has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end




a = param.a;

g = cell(2,1);
g{1} = @(x) simoncelli(x*(2/lmax),a);
g{2} = @(x) real(sqrt(1-(simoncelli(x*(2/lmax),a)).^2));

end


function y = simoncelli(val,a)

y = zeros(size(val));

l1 = a;
l2 = 2*a;

r1ind = val >= 0     &    val < l1;
r2ind = val >= l1    &    val < l2;
r3ind = val >= l2;


y(r1ind) = 1;
y(r2ind) = cos(pi/2*log(val(r2ind)/a)/log(2));
y(r3ind) = 0;


end
