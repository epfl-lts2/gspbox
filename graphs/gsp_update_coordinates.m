function G = gsp_update_coordinates(G, coords)
%GSP_UPDATE_COORDINATES Updates the coordinates of a graph structure
%   Usage: G = gsp_update_coordinates(G, coords);
%
%   Input parameters
%         G      : Graph
%         coords : New coordinates
%
%   Output parameters
%         G      : Output graph
%
%   Update the coordinates of a graph structure
%
%   See also: gsp_compute_coordinates

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

G.coords = coords;
G = rmfield(G,'plotting');
G = gsp_graph_default_plotting_parameters(G);

end