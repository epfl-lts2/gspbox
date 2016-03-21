function G = gsp_update_coordinates(G, coords)
%GSP_UPDATE_COORDINATES Updates the coordinates of a graph structure
%   Usage: G = gsp_update_coordinates(G, coords);
%
%   Input parameters
%         G      : Graph
%         coords : New coordinates
%
%   Output parameters
%         G      : Output graph
%
%   Update the coordinates of a graph structure
%
%   See also: gsp_compute_coordinates
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_update_coordinates.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

G.coords = coords;
G = rmfield(G,'plotting');
G = gsp_graph_default_plotting_parameters(G);

end
