function G = gsp_assign_rand_direction(G,p)
% GSP_ASSIGN_RAND_DIRECTION Assign random direction to p percent edges of the undirected graph G
%
%   Usage: G = gsp_assign_rand_direction(G,p)
%
%   Input parameters:
%       G   : Undirected Graph structure
%       p   : Percentage of edges
%   Output parameters:
%       G   : Directed Graph structure
%
%   Assign random direction to p percent edges of the undirected graph G
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_assign_rand_direction.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

% Author: Francesco Grassi
% Date  : 4 July 2014
% Testing: test_operators

if nargin<2
    p = 1;
end
W = G.W;

[i]=find(triu(W));

k = randperm(length(i),round(p*length(i)));

W(i(k))=0;

G = gsp_update_weights(G,W);

end
