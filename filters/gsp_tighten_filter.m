function [ ftighten ] = gsp_tighten_filter(G, filters )
%GSP_TIGHTEN_FILTER Create a function that tighten a filterbank
%   Usage: ftighten = gsp_tighten_filter( filters );
%
%   Input parameters:
%       G           : Graph or maximum eigenvalue
%       filters     : Filters of the filterbank (cell array)
%   Ouput parameters:
%       ftighten    : Inline function
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using::
%
%       G = gsp_estimate_lmax(G);
%

% Author: Nathanael Perraudin, David Shuman
% Date  : 16 June 2014


if isstruct(G)
    if ~isfield(G,'lmax')
            warning('GSP_TIGHTEN_FILTER: has to compute lmax \n')
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end



[~,B] = gsp_filterbank_bounds([0,lmax],filters);

ftighten = @(x) f_tighten(filters,x,B);

end


function output = f_tighten(filters,x,B)
    output = B;

    for i=1:length(filters)
        output=output-(filters{i}(x)).^2;
    end
    output=real(sqrt(output));
end

