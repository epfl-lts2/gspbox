function [G]=gsp_airfoil()
%GSP_AIRFOIL Initialize the airfoil graph
%   Usage:  G=gsp_airfoil();
%
%   Input parameters:
%         none
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_airfoil()' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.limits) of the airfoil mesh graph. All edge weights are equal
%   to 1.
%
%   Example:
%
%          G = gsp_airfoil();
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   See also: gsp_graph
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_airfoil.php

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


% Author : David I Shuman, Nathanael Perraudin
% Date : 15 March 2014

load airfoil
%G.N=4253;
A = sparse(i_inds,j_inds,1,4253,4253);
G.W = (A + A')/2;
G.coords=[x,y];
G.plotting.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];

G.type = 'airfoil';

G.plotting.vertex_size=30;

G = gsp_graph_default_parameters(G);


end


