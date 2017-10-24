function [ G ] = gsp_cube( param )
%GSP_CUBE Create a graph corresponding to the sampling of an hyper-cube
%   Usage :  G = gsp_cube();
%            G = gsp_cube( param );
%
%   Input parameters:
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_cube( param )' creates a graph from points sampled on a
%   hyper-cube. The dimension of the cube can be passed as a parameter.
%   It can be sampled in a uniform voxel grid or randomly.  
%
%   Additional parameters
%   ---------------------
%
%   * *param.radius*    : float                 the edge length
%   * *param.nb_pts*    : int                   the number of vertices
%   * *param.nb_dim*    : int                   the dimension
%   * *param.sampling*  : ['random']            the variance of the distance kernel
%
%   Example:::
%
%          G = gsp_cube();
%          gsp_plot_graph(G);
%          axis square
%
%   See also: gsp_sphere
%

% Author : Johan Paratte
     if nargin < 1
       param = {};
    end
    
    %Parameters
    if ~isfield(param, 'radius'), param.radius = 1; end
    if ~isfield(param, 'nb_pts'), param.nb_pts = 300; end
    if ~isfield(param, 'nb_dim'), param.nb_dim = 3; end
    if ~isfield(param, 'sampling'), param.sampling = 'random'; end
    
    N = param.nb_pts;
    K = param.nb_dim;
    
    if (K > 3) 
       error('Dimension > 3 not supported yet !'); 
    end
    
    switch param.sampling
        %case 'uniform'
            %Define iteration of each dimension
            
        case 'random'
            % Draw angles randomly
            if (K == 2)
               pts = rand(N, N);
               G = pts;
            end
            
            if (K == 3)
               
               n = floor(N/6);
               
               pts = zeros(n*6, 3);
               pts(1:n, :) = [zeros(n,1), rand(n, 2)];
               pts(n+1:2*n, :) = [ones(n,1), rand(n, 2)];
               
               pts(2*n+1:3*n, :) = [rand(n, 1), zeros(n,1), rand(n, 1)];
               pts(3*n+1:4*n, :) = [rand(n, 1), ones(n,1), rand(n, 1)];
               
               pts(4*n+1:5*n, :) = [rand(n, 2), zeros(n,1)];
               pts(5*n+1:6*n, :) = [rand(n, 2), ones(n,1)];
               
            end

            %TODO Implement http://en.wikipedia.org/wiki/Hypersphere#CITEREFMarsaglia1972
            
        otherwise
            error('Unknown sampling !');
    end
    
    param.type = 'knn';
    param.k = 10;
    G = gsp_nn_graph(pts, param);
    
end

