function gw = gsp_warp_filter(g,w)
%GSP_WARP_FILTER Warp the filterbank g with the filter w
%   Usage: gw = gsp_warp_filter(g,w);
%
%   Input parameters:
%       g   : filterbank
%       w   : warping filter
%
%   Output parameters:
%       gw  : warped filterbank
%
%   The resulting filter *gw* is $gw(x)=w(g(x))$.
%

% Author: Nathanael Perraudin
% Date  : 30 September 2015

Nf = numel(g);

if ~iscell(g)
    g = {g};
end

gw = cell(Nf,1);

for ii = 1:Nf
    gw{ii} = @(x) w(g{ii}(x));
end


end