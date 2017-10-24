function [ bool ] = gsp_test_duality(G, g,h,tol )
%GSP_TEST_DUALITY Test if two filterbanks are dual
%   Usage: bool = gsp_test_duality(G, g,h )
%          bool = gsp_test_duality(G, g,h,tol )
%
%   Input parameters:
%       G       : Graph or arange (min and max value) 
%       g       : filter 1 (or filterbank)
%       h       : filter 2 (or filterbank)
%       tol     : tolerance for the test (default 1e-8)
%
%   Ouput parameters:
%       bool    : boolean 
%
%   This function tests if two filterbanks are dual.
%
%   Example:::
%
%             N = 100;
%             G = gsp_sensor(N);
%             G = gsp_estimate_lmax(G);
%             g = gsp_design_abspline(G,8);
%             gd = gsp_design_can_dual(g);
%             gsp_test_duality(G, g,gd )
%
%   See also: gsp_design_can_dual gsp_test_duality_coefficient

% Author: Nathanael Perraudin
% Date  : 30 December 2014
% Testing: test_dual

if nargin<4
    tol = 1e-8;
end

if isstruct(G)
    if ~isfield(G,'lmax');
        G = gsp_estimate_lmax(G);
        if param.verbose
        warning(['GSP_TEST_DUALITY: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
        end
    end
  arange = [0, G.lmax];
else
  arange = G;
end

N = 100;

x = linspace(arange(1),arange(2),N)';
c1 = gsp_filter_evaluate(g,x);
c2 = gsp_filter_evaluate(h,x);

bool = gsp_test_duality_coefficient(c1,c2,tol);


end

