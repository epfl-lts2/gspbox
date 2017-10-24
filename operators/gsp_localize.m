function [ gt ] = gsp_localize(G, g, n,param)
%GSP_LOCALIZE Localize a kernel g to the node n
%   Usage: gt = gsp_localize(G, g, n);
%
%   Input parameters
%       G   : Graph
%       g   : kernel (or filterbank)
%       n   : Indices of vertex (int)
%       param: Optional parameters
%   Output parameters
%       gt  : translate signal
%
%   This function localize the kernel $g$ onto the node i. If *g*
%   is a cell array, the localization will be done to each filter.
%

% Author: Nathanael Perraudin
% Date  : 28 July 2014

if nargin <4
    param = struct;
end

f = zeros(G.N,numel(n));
for ii = 1:numel(n)
    f(n(ii),ii) = 1;
end
gt = gsp_filter_analysis(G,g,f,param);

end

