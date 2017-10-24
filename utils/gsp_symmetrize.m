function W = gsp_symmetrize(W, type)
%GSP_SYMMETRIZE symmetrize a matrix
%   Usage:  W = gsp_symmetrize(W)
%           W = gsp_symmetrize(W, type)
%
%   Input parameters:
%       W       : square matrix
%       type    : type of symmetrization (default 'full')
%   Output parameters:
%       W       : symmetrized matrix
%
%   The available symmetrization types are:
%   * 'average' : average of $W$ and $W^T$ (default)
%   * 'full'    : copy the missing entries
%   * 'none'    : nothing is done (the matrix might stay unsymmetric!)
%

% Author: Nathanael Perraudin
% Date  : 17 Janvier 2015

if nargin < 2
    type = 'full';
end

switch type
    case 'average'
        W = (W+W.')/2;
    case 'full'
% Old version, slower
%         A = W>0;
%         M = logical(A - (A' & A));
% 
%         W = W + M'.*W';
        W = max(W, W');
    case 'none'
        return
    otherwise
        error('GSP_SYMMETRIZE: Unknown type')
end

end
