function gm = gsp_multiply_filters(g1,g2)
%GSP_MULTIPLY_FILTERS Mutliply to filters
%   Usage: gm = gsp_multiply_filters(g1,g2);
%
%   Input parameters:
%       g1   : filterbank
%       g2   : filterbank
%
%   Output parameters:
%       gm  : multiplied filterbank
%
%   The resulting filter is $gm(x) = g1(x) g2(x)$.
%

% Author: Nathanael Perraudin
% Date  : 30 September 2015

Nf1 = numel(g1);
Nf2 = numel(g2);

if ~iscell(g1)
    g1 = {g1};
end
if ~iscell(g2)
    g2 = {g2};
end

gm = cell(Nf1,Nf2);

for ii = 1:Nf1
    for jj = 1:Nf2
        gm{ii,jj} = @(x) g1{ii}(x).*g2{jj}(x);
    end
end


end