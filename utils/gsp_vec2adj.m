function W = gsp_vec2adj(A,x,type)
%GSP_VEC2ADJ Create the matrix W with sparsity pattern A and entries x
%
%   Usage: W = gsp_vec2adj(A,x)
%              gsp_vec2adj(A,x,type)
%
%   Input parameters:
%       A       : Matrix of sparsity pattern (e.g. binary adjacency matrix)
%       x       : Vector of entries
%       type    : Type of matrix (string) (default 'sym' if A is symmetric)
%                    *  sym: W will be symmetric (x must be of size nnz(A)/2)
%                    * asym: W is asymmetric (x must be of size nnz(A))
%   Output parameters
%       W       : Weighted matrix
%
%   Create the matrix W with sparsity pattern A and entries x%
%

% Author: Francesco Grassi
% Date   : July 2016


if nargin<3
    if issymmetric(full(double(A)))
        type='sym';
    else
        type='asym';
    end
end

if ~nnz(A)==length(x) && ~nnz(A)/2==length(x)
    error('x must have the same size of nnz(A) or nnz(A)/2 if A is symmetric')
end

[N1,N2] = size(A);

switch type
    case 'sym'
        [i,j] = find(triu(A));
        W = sparse(i,j,x,N1,N2)+sparse(i,j,x,N1,N2)';
    case 'asym'
        [i,j] = find(A);
        W = sparse(i,j,x,N1,N2);
    otherwise
        error('Unknow type');
end

    



