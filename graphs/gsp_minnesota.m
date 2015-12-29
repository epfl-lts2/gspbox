function [G] = gsp_minnesota(connect)
%GSP_MINNESOTA  Initialize the Minnesota road network
%   Usage:  G = gsp_minnesota();
%           G = gsp_minnesota(connect);
%
%   Input parameters:
%         connect : change the graph so that it is connected (default 1)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_minnesota()' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.limits) of the Minnesota road network from the MatlabBGL library.
%
%   Remark: if connect is set to 1. We adjust the adjacency matrix so that
%   all edge weights are equal to 1, and the graph is connected.  It is the
%   default!
%
%   To get the orinial disconnected graph, use:
%
%           G = gsp_minnesota(connect);
%
%   Example:
%
%          G = gsp_minnesota();
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   References:
%     D. Gleich. The MatlabBGL Matlab library.
%     http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/index.html.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_minnesota.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

if nargin < 1
   connect = 1; 
end

Q=load('minnesota.mat');
G.N=size(Q.A,1);
G.coords=Q.xy;
G.plotting.limits=[-98,-89,43,50];

if connect
    % Edit adjacency matrix
    A=Q.A;
    % clean minnesota graph
    A=A-diag(diag(A));
    % missing edge needed to connect graph
    A(349,355)=1;
    A(355,349)=1;
    % change a handful of 2 values back to 1
    A(86,88)=1;
    A(88,86)=1;
    A(345,346)=1;
    A(346,345)=1;
    A(1707,1709)=1;
    A(1709,1707)=1;
    A(2289,2290)=1;
    A(2290,2289)=1;
    G.W=sparse(A);
    G.type = 'minnesota';
else
    G.W=sparse(Q.A);
    G.type = 'minnesota-disconnected';
end

G.plotting.vertex_size=30;

G = gsp_graph_default_parameters(G);

end


