function G = gsp_sphere(N, param)
%GSP_SPHERE Create a spherical-shaped graph
%   Usage :  G = gsp_sphere();
%            G = gsp_sphere( param );
%
%   Input parameters:
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_sphere( param )' creates a graph from points sampled on a
%   hyper-sphere. The dimension of the sphere can be passed as a parameter.
%   It can be sampled in a uniform voxel grid or randomly.  
%
%   Additional parameters
%   ---------------------
%
%   * *param.radius*    : float                 the radius of the sphere
%   * *param.nb_pts*    : int                   the number of vertices
%   * *param.nb_dim*    : int                   the dimension
%   * *param.sampling*  : ['random'] the variance of the distance kernel
%
%   Example:::
%
%          G = gsp_sphere();
%          gsp_plot_graph(G);
%          axis square
%
%   See also: gsp_cube
%

% Author : Johan Paratte

if nargin < 1
    N = 300;
end
    

    if nargin < 2
       param = {};
    end
    
    %Parameters
    if ~isfield(param, 'radius'), param.radius = 1; end
    if ~isfield(param, 'nb_dim'), param.nb_dim = 3; end
    if ~isfield(param, 'sampling'), param.sampling = 'random'; end

    K = param.nb_dim;
    
    switch param.sampling
        case 'uniform'
            % Recursive zonal equal area partition EQ(d, N)
            pts = eq_point_set(K - 1, N)';
            
        case 'random'
            % Draw angles randomly
%             angles = rand(N, K-1);
%             if (K > 2)
%                 angles(:,1:K-2) = angles(:,1:K-2)*2*pi;
%                 angles(:,K-1) = angles(:,K-1)*pi;
%             else 
%                 angles(:,K-1) = angles(:,K-1)*2*pi;
%             end
%             angles(:,K-1) = angles(:,K-1)*pi;
%             pts = ones(N, K);
%             % x_1 = r
%             pts(:,1) = pts(:,1) .* param.radius;
%             for k = 2:K-1
%                pts(:,k) = pts(:,k-1) .* sin(angles(:,k-1)); 
%             end
%             pts(:,K) = pts(:,K-1);
%             pts(:,1:K-1) = pts(:,1:K-1) .* cos(angles);
%             pts(:,K) = pts(:,K) .* sin(angles(:,K-1));
            
            pts = normrnd(0, 1, N, K);
            for ii = 1:N
               pts(ii,:) = pts(ii,:) ./ norm(pts(ii,:));
            end

            %TODO Implement http://en.wikipedia.org/wiki/Hypersphere#CITEREFMarsaglia1972
            
        otherwise
            error('Unknown sampling !');
    end
    
    param.type = 'knn';
    param.k = 10;
    G = gsp_nn_graph(pts, param);

end

