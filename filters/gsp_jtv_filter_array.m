function [h,filtertype] = gsp_jtv_filter_array(G,g,filtertype)
%GSP_JTV_FILTER_ARRAY Convert ts/js filters to -array filters
%   Usage: [h,filtertype] = gsp_jtv_filter_array(G,g, filtertype)
%
%   Input parameters:
%         G          : Graph
%         g          : Cell array of time-vertex filters
%         filtertype : Filter domain (ts,js)
%   Output parameters:
%         h          : Cell array of graph filterbank
%         filtertype : Filter domain (ts-array,js-array)
%
%   Convert ts/js filters to -array filters
%

% Author :  Francesco Grassi, Nathanael Perraudin
% Date : September 2016

if ~iscell(g)
   g = {g};
end

if ~gsp_check_filtertype(filtertype,{'ts','js'})
    error('Invalid filtertype.');
end

T  = G.jtv.T;
Nf = numel(g);

switch filtertype
    case 'ts'
        v = gsp_jtv_ta(G);
    case 'js'
        v = gsp_jtv_fa(G);
end

h  = cell(Nf,T);
for n=1:Nf
    for ii = 1:T
        h{n,ii} = @(x) g{n}(x,v(ii));
    end
end

filtertype = [filtertype '-array'];

end