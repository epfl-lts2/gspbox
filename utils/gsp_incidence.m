function G = gsp_incidence(G,inc_type)
%GSP_INCIDENCE Compute an incidence matrix
%
%   Usage: G = gsp_incidence(G)
%              gsp_incidence(G,inc_type)
%
%   Input parameters:
%       G        : Graph structure
%       inc_type : Type of incidence matrix: 'weighted' or 'binary'
%   Output parameters:
%       G   : Graph structure
%
%   'gsp_incidence' compute the incidence matrix B from the gradient and add it to graph structure
%

% Author: Francesco Grassi
% Date  : July 2016

if nargin < 2
    inc_type = 'weighted';
end


if ~isfield(G,'Diff'); G = gsp_adj2vec(G); end;

switch inc_type
    case 'weighted'
        B = abs(G.Diff);
    case 'binary'
        B = abs(G.Diff)>0;
    otherwise
        error('Unknown incidence matrix.')
end

G.B = B;

end