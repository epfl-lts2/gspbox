function Gn = gsp_copy_graph_attributes(G,type,Gn)
%GSP_COPY_GRAPH_ATTRIBUTES copy the parameter of the graph
%   Usage: Gn = gsp_copy_graph_attributes(G);
%
%   Input arguments:
%       G       : Graph structure
%       type    : flag to select what to copy (default 1)
%       Gn      : Graph structure (optional)
%
%   Output arguments:
%       Gn      : Partial graph structure
%
%   This function copy optional argument of a graph but not the graph
%   itself. If a graph is given as a third argument, then it will copy
%   everything in this graph. Otherwise, the function will create a new
%   empty graph.
%
%
%   The flag type can take the following value:
%
%    0: copy only the the graph structure
%        lap_type
%        plotting (but not plotting.limits)
%
%    1: copy additionaly the field
%        plotting.limits
%        coords
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_copy_graph_attributes.php

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
% Date  : 24 July 2014


if nargin<2
    type = 1;
end

if nargin<3
    Gn = struct;
end

if isfield(G, 'lap_type'), Gn.lap_type = G.lap_type; end
if isfield(G, 'plotting'), Gn.plotting = G.plotting; end



if type
    if isfield(G, 'coords'), Gn.coords = G.coords; end
else
    if isfield(Gn.plotting,'limits')
        Gn.plotting = rmfield(Gn.plotting,'limits');
    end
end


end
