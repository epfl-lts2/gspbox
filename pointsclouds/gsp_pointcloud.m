function [P, info] = gsp_pointcloud( name, max_dim)
%GSP_POINTCLOUD Load models and return the points
%   Usage:  P = gsp_pointcloud(name)
%           P = gsp_pointcloud(name, max_dim)
%
%   Input parameters:
%           name        : the name of the point cloud to load ('airfoil', 'two_moons', 'bunny')
%           max_dim     : the maximum dimensionality of the points (only valid for two_moons)
%
%   Output parameters:
%           P           : set of points in a NxD with N the number of points and D the dimensionality of the pointcloud
%           info        : optional additional information
%
%   'gsp_pointcloud( name, max_dim)' load pointcloud data and format it in
%   a unified way as a set of points with each dimension in a different
%   column
%
%   Note that the bunny is the model from the Stanford Computer Graphics
%   Laboratory see references. 
%
%   See also: gsp_nn_graph
%
%   References:
%     G. Turk and M. Levoy. Zippered polygon meshes from range images. In
%     Proceedings of the 21st annual conference on Computer graphics and
%     interactive techniques, pages 311--318. ACM, 1994.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/gsp_pointcloud.php

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

% Author: Johan Paratte
% Date: 7 August 2014
% 

    if nargin < 2
    % Select different defaults depending on the point cloud
        max_dim = 2;
    end

     switch name
        % Load the airfoil and fill the points
        case 'airfoil'
            data = load('airfoil.mat');
            info.i_inds = data.i_inds;
            info.j_inds = data.j_inds;
            P = [data.x, data.y];
        % Load the two moons and fill the points
        case 'two_moons'
            data = load('two_moons.mat');
            if (max_dim == -1) 
                max_dim = 2;
            end
            P = data.features(1:max_dim, :)';
        % Load the bunny model and fill the points
        case 'bunny'
             data = load('bunny.mat');
             P = data.bunny;
        otherwise
            error('Unknown pointcloud !');
    end

end


