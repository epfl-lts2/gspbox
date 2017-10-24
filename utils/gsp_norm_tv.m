function y = gsp_norm_tv(G,x)
%GSP_NORM_TV TV norm on graph
%   Usage:  y = gsp_norm_tv(G,x);
%
%   Input parameters:
%         G     : Graph structure
%         x     : Signal on graph
%   Output parameters:
%         y     : Norm
%
%   Compute the TV norm of a signal on a graph
%
%   See also: gsp_prox_tv

% Author: Nathanael Perraudin
% Date:   25 March 2014
% Testing: test_gsp_prox

if ~isfield(G,'v_in')
    G = gsp_adj2vec(G);
    warning(['GSP_PROX_TV: To be more efficient you should run: ',...
        'G = gsp_adj2vec(G); before using this proximal operator.']);
end


y = sum(abs(gsp_grad(G,x)));

end
