function C = gsp_strong(A,B)
%GSP_STRONG Strong product between matrices
%   Usage:  C = gsp_strong(A,B)
%
%   Input parameters:
%         A          : matrix size N1*M1
%         B          : matrix size N2*M2
%   Output parameters:
%         C          : matrix size (N1xN2) * (M1xM2)
%
%   'gsp_strong' computes the strong product between two matrices A
%   and B. The cartesian product between two matrices is
%
%   .. C = strong(A,B) = kron(A,B) + cartesian(A,B) 
%
%   .. math:: C = A \boxtimes B = A \otimes B +  A \otimes I_{N_2} + I_{N_1} \otimes B
%
%   See also: gsp_graph_product

% Author : Francesco Grassi
% Date   : July 2016

C = kron(A,B) + gsp_cartesian(A,B);


end