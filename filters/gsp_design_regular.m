function [ g ] = gsp_design_regular(G, param)
%GSP_DESIGN_REGULAR Create a regular filterbank
%   Usage: g = gsp_design_regular( G );
%          g = gsp_design_regular( G, param );
%   
%   Inputs parameters:
%       G       : Graph structure or lmax
%       param   : Structure of optional parameters
%
%   Outputs parameters:
%       g       : filterbank
%
%   This function creates a parseval filterbank of $2$ filters. The low-pass
%   filter is defined by a function $f_l(x)$ between $0$ and $2$. For
%   $d = 0$.
%
%   ..    f_l(x) = sin(pi/4*x)
%
%   .. math:: f_{l}= \sin\left( \frac{\pi}{4} x \right)
%
%   For $d = 1$ 
%
%   ..    f_l(x) = sin( pi/4 * (1+sin(pi/2*(x-1))) )
%
%   .. math:: f_{l}= \sin\left( \frac{\pi}{4} \left( 1+ \sin\left(\frac{\pi}{2}(x-1)\right) \right) \right)
%
%   For $d = 2$ 
%
%   ..    f_l(x) = sin( pi/4 * ( 1 + sin( pi/2 * sin(pi/2*(x-1) ) ) )
%
%   .. math:: f_{l}= \sin\left( \frac{\pi}{4} \left( 1+ \sin\left(\frac{\pi}{2} \sin\left(\frac{\pi}{2}(x-1)\right)\right) \right) \right)
%
%   And so on for the other degrees $d$.
%
%   The high pass filter is adapted to obtain a tight frame.
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
%         g = gsp_design_regular(G);   
%         gsp_plot_filter(G,g);  
%         [A,B] = gsp_filterbank_bounds(G,g)
%
%   *param* is an optional structure containing the following fields
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.d*: Degree. See equation for mor informations. (default 3)
%

% Author: Nathanael Perraudin, David Shuman
% Date  : 21 June 2014
% Testing: test_filter


if nargin < 2
    param = struct;
end


if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'d'), param.d = 3; end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_REGULAR has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end




d = param.d;

g = cell(2,1);
g{1} = @(x) regular(x*(2/lmax),d);
g{2} = @(x) real(sqrt(1-(regular(x*(2/lmax),d)).^2));

end


function y = regular(val,d)


if d==0
    y = sin(pi/4*val);
else
    output = sin(pi*(val-1)/2);
    for k=2:d
        output = sin(pi*output/2);
    end
    y = sin(pi/4*(1+output));
end


end
