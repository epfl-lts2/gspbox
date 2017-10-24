function [G]=gsp_logo()
%GSP_LOGO Initialize a graph with the GSP logo
%   Usage:  G = gsp_logo();
%
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_logo()' initializes a graph structure containing the GSP logo
%
%   Example:
%
%          G = gsp_logo();
%          gsp_plot_graph(G);
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_logo.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

% Author : Johan Paratte
% Test: 

load logogsp

G.W = W;
G.coords = coords;
G.info.idx_g = idx_g;
G.info.idx_s = idx_s;
G.info.idx_p = idx_p;
%G.plotting.vertex_color = [200 136/255.0 204/255.0];
%G.plotting.edge_color = [0 136/255.0 204/255.0];
G.plotting.vertex_size = 20;
G.limits = [0 640 -400 0];
G = gsp_graph_default_parameters(G);

end



