function n = gsp_norm_l2_filterbank(G, W, x, param)
%GSP_NORM_L2_FILTERBANK Compute the l2 norm of the analysis coefficients
%   Usage: n = gsp_norm_l2_filterbank(G, W, x);
%
%   Input parameters:
%       G   : Graph structure
%       W   : Filterbank (cell array of functions)
%       x   : coefficients
%       param: structure of optional paramter
%   Output parameters:
%       n   : L2 norm
%
%   `gsp_norm_l2_filterbank(G, W, x, param)` computes:
%
%   .. n = || A W^* x -y ||_2^2
%
%   .. math::  n =  \| A W^* x - y\|_2^2 
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.A* : Forward operator (default: Id).
%
%   * *param.y* : measurements (default: 0).
%

% Author: Nathanael Perraudin
% Date: 26 March 2014

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
warning(['GSP_NORM_L1_FILTERBANK: To be more efficient you should run: ',...
    'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

if nargin < 4, param=struct; end

if ~isfield(param, 'y'), param.nu = 0; end
if ~isfield(param, 'A'), param.A = @(x) x; end

n = norm(param.A(gsp_filter_synthesis(G,W,x))-param.y,'fro')^2;

end
