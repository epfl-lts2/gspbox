function [ Gout ] = gsp_separate_graph( G )
%GSP_SEPARATE_GRAPH Separate the graph G into disconnected components
%   Usage: [ Gout ] = gsp_separate_graph( G )
%
%   Input parameters:
%       G   : Original graph
%
%   Output parameters:
%       Gout: Cell array of graphs
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_separate_graph.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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
% Date : Mai 2014

% TODO: DO this operation better
G = gsp_compute_fourier_basis(G);

NG = sum(G.E < 1e-10);

error('Not implemented');

Gout = 0;


end


