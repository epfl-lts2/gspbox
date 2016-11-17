function [I, J] = ind2sub4up(IND)
%IND2SUB4UP Subscripts from linear index for upper triangular matrix (only
%elements above diagonal)
%   IND2SUB4UP determines the equivalent subscript values corresponding to
%   a given single index into a 2D upper triangular matrix, excluded all
%   elements over the diagonal.
%
%   [I, J] = IND2SUB4UP(IND) returns vectors I and J containing equivalent
%   row and column subscripts corresponding to the index vector IND.
%
%   Let ind be a vector of indexes for entries of some upper triangular
%   matrix. The entries are selected vertically so that:
%
%       ind = 1                is associated to entry      (1, 2)
%       ind = 2                is associated to entry      (1, 3)
%       ind = 3                is associated to entry      (2, 3)
%       ind = 4                is associated to entry      (1, 4)
%       ...
%       ind = N * (N - 1) / 2  is associated to entry      (N - 1, N)
%
% % ======================================================================
%
%    EXAMPLE
%
%       % Note that if
%             A = rand(10);
%       % and
%             b = A(find(triu(A, 1)));
%       % then, given indices
%             IND = [1:45];
%       % for vector b, these are equivalent to subscripts
%             [I, J] = ind2sub4up(IND);
%       % for matrix A. In fact:
%             all(A(sub2ind(size(A), I, J)) == b(IND))
%
%       %    ans =
%       %           1
%
%       % This is obtained without even knowing about size(A)
%
% % ======================================================================
%
%   See also SUB2IND, IND2SUB, FIND.
%
% % ======================================================================
%
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%                                                                                               %
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 1 May 2010                                                       %
%                                                                                               %
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%
% % ======================================================================
%

% Check input
ctrl1 = isnumeric(IND) & isreal(IND);
if ctrl1
  IND = ceil(IND(:));
  ctrl2 = ~any(isnan(IND)) & ~any(isinf(IND)) & all(IND > 0);
  if ~ctrl2
  error('Check indexes: they need to be positive integers!')
  end
else
  error('Check indexes: they need to be positive integers!')
end

J = round(floor(-.5 + .5 * sqrt(1 + 8 * (IND - 1))) + 2);
I = round(J .* (3 - J) / 2 + IND - 1);
