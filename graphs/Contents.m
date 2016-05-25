% GSPBOX - Graphs
%
%  Specific graphs
%    gsp_swiss_roll              -  Create swiss roll graph
%    gsp_david_sensor_network    -  Create the sensor newtwork from david
%    gsp_ring                    -  Create the ring graph
%    gsp_path                    -  Create the path graph
%    gsp_airfoil                 -  Create the airfoil graph
%    gsp_comet                   -  Create the comet graph
%    gsp_erdos_renyi             -  Create a erdos renyi graph
%    gsp_minnesota               -  Create Minnesota road graph
%    gsp_low_stretch_tree        -  Create a low stretch tree graph
%    gsp_sensor                  -  Create a random sensor graph
%    gsp_random_regular          -  Create a random regular graph
%    gsp_random_ring             -  Create a random ring graph
%    gsp_full_connected          -  Create a fully connected graph
%    gsp_nn_graph                -  Create a nearest neighbors graph
%    gsp_rmse_mv_graph           -  Create a nearest neighbors graph with missing values
%    gsp_sphere                  -  Create a spherical-shaped graph
%    gsp_cube                    -  Create a cubical-shaped graph
%    gsp_2dgrid                  -  Create a 2d-grid graph
%    gsp_torus                   -  Create a torus graph
%    gsp_logo                    -  Create a GSP logo graph
%    gsp_community               -  Create a community graph
%    gsp_bunny                   -  Create a bunny graph
%    gsp_spiral                  -  Create a spiral graph
%    gsp_stochastic_block_graph  -  Create a graph with the stochastic block model
%
%  Hypergraphs
%    gsp_nn_hypergraph           -  Create an hyper nearest neighbor graph
%    gsp_hypergraph              -  Create an hypergraph
%
%  Utils
%    gsp_graph_default_parameters-  Initialise all parameters for a graph
%    gsp_graph_default_plotting_parameters-  Initialise all plotting parameters for a graph
%    gsp_graph                   -  Create a graph from a weight matrix
%    gsp_update_weights          -  Update the weights of a graph
%    gsp_update_coordinates      -  Update the coordinate of a graph
%    gsp_components              -  Cuts non connected graph into several connected ones
%    gsp_subgraph                -  Create a subgraph
%
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/Contents.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

