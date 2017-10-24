function [A, B] = gsp_filterbank_bounds(G,W,param)
%GSP_FILTERBANK_BOUNDS Compute approximated frame bounds for a filterbank
%   Usage: [A, B] = gsp_filterbank_bounds(G,W);
%          [A, B] = gsp_filterbank_bounds(G,W,param);
%          [A, B] = gsp_filterbank_bounds([xmin, xmax],W);
%          [A, B] = gsp_filterbank_bounds([xmin, xmax],W,param);
%
%   Input parameters
%       G   : Graph structure or interval to compute the bound
%       W   : Filterbank (cell array of inline function)
%       param: optional parameter
%   Output parameters
%       A   : Filterbank lower bound
%       B   : Filterbank Upper bound
%
%   *param* is a Matlab structure containing the following fields:
%
%   * *param.N* : Number of points for the line search default (default 999)
%   * *param.use_eigenvalues* : Use eigenvalues if possible (default 1). To
%     be used, the eigenvalues have to be computed first using
%     |gsp_compute_fourier_basis|.
%

% Author: Nathanael Perraudin
% Date : 26 March 2014

if nargin < 3
    param = struct;
end

if ~isfield(param,'N'), param.N = 999; end
if ~isfield(param,'use_eigenvalues'),  param.use_eigenvalues = 1; end

if iscell(G)
    NG = numel(G);
    A = cell(NG,1);
    B = cell(NG,1);
    for ii = 1:NG
        [A{ii}, B{ii}] = gsp_filterbank_bounds(G{ii},W{ii},param);
    end
    return
end

if isstruct(G)
    if ~isfield(G,'lmax')
        G = gsp_estimate_lmax(G);
    warning(['GSP_FILTERBANK_BOUNDS: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
    end
    xmax = G.lmax;
    xmin = 0;
else
    xmin = G(1);
    xmax = G(2);
end
    

if param.use_eigenvalues && isstruct(G) && isfield(G,'e')
    lambda = G.e;
else
    lambda = linspace(xmin,xmax , param.N);
end

 
Nf = numel(W);


sum_filters = sum(abs(gsp_filter_evaluate(W,lambda).^2),2);
% for ii=1:Nf
%     sum_filters = sum_filters + (W{ii}(lambda)).^2;
% end

A = min(sum_filters);
B = max(sum_filters);
  
  
  
end
