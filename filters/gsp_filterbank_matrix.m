function [ F ] = gsp_filterbank_matrix(G,g,param )
%GSP_FILTERBANK_MATRIX Create the matrix of the filterbank frame
%   Usage:  F = gsp_filterbank_matrix(G, g param );
%           F = gsp_filterbank_matrix(G,g);
%
%   Input parameters:
%         G     : Graph
%         g     : Filters
%         param : Structure of optional parameter
%   Output parameters:
%         F     : Frame
%
%   This function creates the matrix associated to the filterbank g. The
%   size of the matrix is MN x N, where M is the number of filters.
%
%   *param* a Matlab structure containing the following fields:
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%



%   AUTHOR : Nathanael Perraudin
%   TESTING: test_filter

% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end

if param.verbose && G.N>2000
    warning('Create a big matrix, you can use other methods.');
end

Nf = length(g);

Ft = gsp_filter_analysis(G,g,eye(G.N));

F = zeros(size(Ft'));
for ii = 1:Nf
    F(:,G.N*(ii-1)+(1:G.N)) = Ft(G.N*(ii-1)+(1:G.N),:);
end

end

