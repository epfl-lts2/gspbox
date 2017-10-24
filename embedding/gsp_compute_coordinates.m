function coords = gsp_compute_coordinates(G, dim, method, param)
%GSP_COMPUTE_COORDINATES calculate coordinates
%   Usage: coords = gsp_compute_coordinates( G , dim );
%          coords = gsp_compute_coordinates( G , dim, method );
%          coords = gsp_compute_coordinates( G , dim, method, param );
%
%   Input parameters
%     G       : Graph
%     dim     : Number of lower dimentions 
%     method  : Method to be used
%     param   : Structure of optional paramter 
%
%   Output parameters
%     coords   : New coordinates of graph that lives in a lower dimention
%    
%   This function computes coordinates from the weighted adjacency matrix. 
%   The available methods are:
%
%   * 'laplacian_eigenmaps' : Laplacian Eigenmaps. For further information,
%     see the help of |gsp_laplacian_eigenmaps|  (default)
%   * 'lle'                 :  LLE (Loccaly Linear Embedding). For further
%     information, see the help of |gsp_lle|
%   * 'isomap'              :  Isomap. For further information, see the
%     help of |gsp_isomap| 
%
%   The coordinates can be insterted inside the graph structure with the
%   function |gsp_update_coordinates|.
%
%   *param* is a structure forwarding optiona paramter to the different
%   functions.
%
%   * *param.tol* : Tolerance for the spectral gap (default 1e-6).   
%     For further parameters, check the selected method function's help.
%
%   See also: gsp_update_coordinates gsp_isomap gsp_laplacian_eigenmaps gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

if nargin<4
    param     = struct;
end

if nargin<3
    method = 'laplacian_eigenmaps';
end

if nargin<2
    dim = 3;
end

if ~isfield(param,'tol'), param.tol = 1e-6; end

if iscell(G)
    coords = cell(size(G));
    for ii = 1: numel(G)
        coords = gsp_compute_coordinates(G{ii}, dim, method, param);
    end 
end


if gsp_check_connectivity(G)==0
    error('Graph is not connected')
    warning('Use gsp_components to extract individual graph components')
    
    [~, Gz, ~] = gsp_components(G);
%     coords = cell(size(Gz));
%     for ii=1:length(Gz)
%         switch lower(method)
%             case 'laplacian_eigenmaps'
%                 coords{ii} = gsp_laplacian_eigenmaps(Gz{ii}, dim, param);
%             case 'lle'
%                 coords{ii} = gsp_lle(Gz{ii}, dim, param);
%             case 'isomap'
%                 coords{ii} = gsp_isomap(Gz{ii}, dim);
%             otherwise
%                 warning('Please specify param.method')
%         end
%    end
    
    coords = gsp_compute_coordinates(G, dim, method, param);

else
    switch lower(method)
        case 'laplacian_eigenmaps'
            coords = gsp_laplacian_eigenmaps(G, dim, param);
        case 'lle'
            coords = gsp_lle(G, dim, param);
        case 'isomap'
            coords = gsp_isomap(G, dim);
        otherwise
            warning('Please specify param.method')
    end
end

end

