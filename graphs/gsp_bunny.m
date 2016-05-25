function [ G ] = gsp_bunny()
%GSP_BUNNY Create a graph of the stanford bunny
%   Usage :  G = gsp_bunny();
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_bunny()' creates a graph from the pointcloud of the Stanford Bunny model. 
%
%   Example:
%
%           G = gsp_bunny();
%           gsp_plot_graph(G);
%
%   References:
%     G. Turk and M. Levoy. Zippered polygon meshes from range images. In
%     Proceedings of the 21st annual conference on Computer graphics and
%     interactive techniques, pages 311--318. ACM, 1994.
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_bunny.php

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

    %Load the point cloud
    P = gsp_pointcloud('bunny');
    
    %Create the graph from the point cloud using an epsilon-neighborhood
    %connectivity
    param.type = 'radius';
    param.rescale = 1;
    param.center = 1;
    param.epsilon = 0.2;
    
    %Compute it
    G = gsp_nn_graph(P, param);
    %Reduce vertex size for plotting
    G.plotting.vertex_size = 10;
end


