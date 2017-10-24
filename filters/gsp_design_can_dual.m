function [ gd ] = gsp_design_can_dual( g ,tol)
%GSP_DESIGN_CAN_DUAL This function return the canonical dual filters of g
%   Usage:  gd = gsp_design_can_dual( g );
%
%   Inputs parameters:
%       g       : cell array of filters
%       tol     : tolerance for the pseudo-inverse
%
%   Ouputs parameters:
%       g       : cell array of filters
%
%   This function returns the canonical dual filterbank g. Note that it
%   might not be the be the optimal solution in term of computation.
%
%   Example:::
%
%             N = 100;
%             G = gsp_sensor(N);
%             G = gsp_compute_fourier_basis(G);
%             g = gsp_design_abspline(G,8);
%             gd = gsp_design_can_dual(g);
%             paramplot.show_sum = 0;
%             figure(1)
%             gsp_plot_filter(G,g,paramplot);
%             title('Original filters')
%             figure(2)
%             gsp_plot_filter(G,gd,paramplot);
%             title('Canonical dual filters');
% 
%             x = rand(N,1);
%             param.method = 'exact';
%             coeff = gsp_filter_analysis(G,g,x,param);
%             xs = gsp_filter_synthesis(G,gd,coeff,param);
%             norm(xs-x)
%
%   See also: gsp_evaluate_can_dual

% Author: Nathanael Perraudin
% Date  : 30 December 2014
% Testing: test_dual

if nargin<2
    tol = 1e-8;
end

Nf = length(g);
gd = cell(Nf,1);

for ii = 1:Nf
    gd{ii} = @(x) can_dual(g,ii,x,tol);
end
    
end


function ret = can_dual(g,n,x,tol)
    [N1, N2] = size(x);
    x = x(:);
    sol = gsp_evaluate_can_dual( g,x,tol );
    ret = sol(:,n);
    ret = reshape(ret,N1,N2);
end
