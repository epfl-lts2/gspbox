function d = gsp_hop_distanz(G,i,j)
%GSP_HOP_DISTANZ Compute the hop distance between two node
%   Usage:  d = gsp_hop_distanz(G,i,j);
%
%   Input parameters:
%       G   : Graph
%       i   : node
%       j   : node
%   Output parameters:
%       d   : hop distanz
%
%   This code computes the hop distance between node i and node j. It uses
%   a naive greedy algorithm and has to be improved.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_hop_distanz.php

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

% Author: Nathanael Perraudin
% Date  : 15 septembre 2015
% Testing: test_gsp_hope_distanz

M = double(logical(G.W));
s = zeros(G.N,1);
s(i) = 1;
s = logical(s);
d = 0;
while s(j)==0 && d<=G.N+1
    d = d+1;
    s = logical(M*double(s)) +s;
end

if d == G.N+1
    d = inf;
end


end
