% sgwt_check_connected : Check connectedness of graph
%
% function r=sgwt_check_connected(A)
% 
% returns 1 if graph is connected, 0 otherwise
% Uses boost graph library breadth first search
%
% Inputs : 
% A - adjacency matrix
%
% Outputs :
% r - result
%

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

function r=sgwt_check_connected(A)
  d=bfs(A,1);
  r=~any(d==-1);
  
