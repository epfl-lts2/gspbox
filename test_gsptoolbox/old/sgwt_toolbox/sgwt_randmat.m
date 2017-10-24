% sgwt_randmat : Compute random  (Erdos-Renyi model) graph
%
% function A=sgwt_randmat(N,thresh)
%
% Inputs : 
% N - number of vertices
% thresh - probability of connection of each edge
%
% Outputs :
% A - adjacency matrix

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function [A]=sgwt_randmat(N,thresh)
  assert(thresh<=1 && thresh>=0);
  A=rand(N)>1-thresh;
  B=triu(A);
  A=B+B';
  for i=1:size(A,1)
    A(i,i)=0;
  end
  A=sparse(A);
