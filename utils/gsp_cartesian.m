function C = gsp_cartesian(A,B)
%GSP_CARTESIAN Cartesian product between vectors or matrices
%   Usage:  C = gsp_cartesian(A,B)
%
%   Input parameters:
%         A          : vector size N1 or matrix size N1 x M1
%         B          : vector size N2 or matrix size N2 x M2
%   Output parameters:
%         C          : vector size N1*N2 or matrix size (N1*N2) x (M1*M2)
%
%   If A and B are vectors 'gsp_cartesian' computes the following product:
%
%   .. C = cartesian(A,B) = kron(A,ones(size(B))) + kron(ones(size(A),B)
%
%   .. math:: C = A \times B = A \otimes 1_{N_2} +  1_{N_1} \otimes B
%
%   If A and B are matrices 'gsp_cartesian' computes the following product:
%
%   .. C = cartesian(A,B) = kron(A,eye(size(B))) + kron(eye(size(A),B)
%
%   .. math:: C = A \times B = A \otimes I_{N_2} + I_{N_1} \otimes B
%
%   See also: gsp_graph_product 

% Author : Francesco Grassi
% Date   : July 2016

%TO DO: extend to more than 2 factors

if isvector(A) && isvector(B)
    C = kron(A,ones(size(B))) + kron(ones(size(A)),B);
elseif ismatrix(A) && ismatrix(B)
    C = kron(A,eye(size(B))) + kron(eye(size(A)),B);
else
    error('GSP_CARTESIAN: Cartesian product can be computed only between matrices or between vectors')
end

end