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
%    'laplacian_eigenmaps' : Laplacian Eigenmaps. For further information,
%     see the help of GSP_LAPLACIAN_EIGENMAPS  (default)
%    'lle'                 :  LLE (Loccaly Linear Embedding). For further
%     information, see the help of GSP_LLE
%    'isomap'              :  Isomap. For further information, see the
%     help of GSP_ISOMAP 
%
%   The coordinates can be insterted inside the graph structure with the
%   function GSP_UPDATE_COORDINATES.
%
%   param is a structure forwarding optiona paramter to the different
%   functions.
%
%    param.tol : Tolerance for the spectral gap (default 1e-6).   
%     For further parameters, check the selected method function's help.
%
%   See also: gsp_update_coordinates gsp_isomap gsp_laplacian_eigenmaps gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/embedding/gsp_compute_coordinates.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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


