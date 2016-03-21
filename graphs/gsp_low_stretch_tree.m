function [G]=gsp_low_stretch_tree(k)
%GSP_LOW_STRETCH_TREE  Initialize a low stretch tree
%   Usage:  G = gsp_low_stretch_tree(k);
%           G = gsp_low_stretch_tree();
%
%   Input parameters:
%         k     : 2^k points on each side of the grid of vertices. (default 6)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_create_low_stretch_tree(k)' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), the plotting coordinate limits 
%   (G.limits), and the root of a low stretch tree on a grid of points.
%   There are 2^k points on each side of the grid, and therefore 2^{2k} 
%   total vertices. The edge weights are all equal to 1.
%
%   Example:
%
%          G = gsp_low_stretch_tree(3);
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   See also: gsp_graph
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_low_stretch_tree.php

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


% Author : David I Shuman, Nathanael Perraudin

if nargin < 1
   k = 6; 
end

ii=[2 3 1 1 4 3];
jj=[1 1 2 3 3 4];

XCoords=[1, 2, 1, 2];
YCoords=[1, 1, 2, 2];

for p=2:k
    % create the indicies of the weight matrice
    ii_new=[ii,ii+4^(p-1),ii+2*4^(p-1),ii+3*4^(p-1)];
    ii_new_middle=[4^(p-1),4^(p-1),4^(p-1)+(4^p+2)/3,5/3*4^(p-1)+1/3,4^(p-1)+(4^p+2)/3,3*4^(p-1)+1];
    ii=[ii_new,ii_new_middle];

    jj_new=[jj,jj+4^(p-1),jj+2*4^(p-1),jj+3*4^(p-1)]; 
    jj_new_middle=[5/3*4^(p-1)+1/3,4^(p-1)+(4^p+2)/3,3*4^(p-1)+1,4^(p-1),4^(p-1),4^(p-1)+(4^p+2)/3];
    jj=[jj_new,jj_new_middle];
    
    % Create Coords
    YCoords=repmat(YCoords,1,2);
    YCoords_new=[YCoords,YCoords+2^(p-1)];
    YCoords=YCoords_new;

    XCoords_new=[XCoords,XCoords+2^(p-1)];
    XCoords=repmat(XCoords_new,1,2);
end

G.W=sparse(ii,jj,ones(size(ii)));
G.coords=[XCoords',YCoords'];

G.limits=[0,2^k+1,0,2^k+1];
G.N=(2^k)^2;
G.root=4^(k-1);

G.plotting.edge_width=1.25;
G.plotting.vertex_size=75;

G.type = 'low strech tree';

G = gsp_graph_default_parameters(G);


end

