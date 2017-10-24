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
%   References: turk1994zippered
%

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

