function [G]=gsp_david_sensor_network(N)
%GSP_DAVID_SENSOR_NETWORK  Initialize a sensor network
%   Usage:  G=gsp_david_sensor_network(N);
%
%   Input parameters:
%         N     : Number of vertices (default 64)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_david_sensor_network(N)' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.limits) of a random sensor network with N vertices. The 
%   sensors are placed randomly in the unit square, and edges are placed 
%   between any sensors within a fixed radius of each other. The edge 
%   weights are assigned via a thresholded Gaussian kernel. The sensor 
%   network will be connected for N=500 or N=64. 
%
%   Warning: this graph is not necessarly connected...
%
%   Example:
%
%          G = gsp_david_sensor_network(64);
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_david_sensor_network.php

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

% Author : David I Shuman, Nathanael Perraudin
% Test: test_graphs

if nargin<1
    N = 64;
end

% TODO: To be changed
    randn('seed', 18); 
    rand('seed', 18);

G.N=N;

if N==64
    load('david64.mat');
    G.W = W;
    G.N = N;
    G.coords = coords;
elseif N==500
    load('david500.mat');
    G.W = W;
    G.N = N;
    G.coords = coords;
else
    % Generate sensor locations
    Xcoords = rand(N,1);
    Ycoords = rand(N,1);
    G.coords = [Xcoords,Ycoords];

    % Create weighted adjancency matrix
    target_dist_cutoff = -.125*N/436.075+.2183;
    T = .6; 
    s = sqrt(-target_dist_cutoff^2/(2*log(T)));
    d = gsp_distanz(G.coords'); 
    G.W = exp(-d.^2/(2*s^2)); 
    G.W(G.W<T) = 0; % Thresholding to have sparse matrix
end

G.plotting.limits = [0,1,0,1];
G.W=G.W-diag(diag(G.W));
G.W=sparse(G.W);

G = gsp_graph_default_parameters(G);

end

