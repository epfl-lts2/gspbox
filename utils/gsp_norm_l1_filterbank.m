function n = gsp_norm_l1_filterbank(G, W, x)
%GSP_NORM_L1_FILTERBANK Compute the l2 norm of the analysis coefficients

%   Usage: n = gsp_norm_l1_filterbank(G, W, x);
%
%   Input parameters:
%       G   : Graph structure
%       W   : Filterbank (cell array of functions)
%       x   : coefficients
%       param: structure of optional paramter
%   Output parameters:
%       n   : L1 norm
%
%   `gsp_norm_l1_filterbank(G, W, x)` computes:
%
%   .. n = || W^* x ||_1
%
%   .. math::  n =  \|  W^* x \|_1
%
%

% Author: Nathanael Perraudin
% Date: 26 March 2014

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
warning(['GSP_NORM_L1_FILTERBANK: To be more efficient you should run: ',...
    'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

n = sum(sum(abs(gsp_filter_analysis(G,W,x))));

end

 
