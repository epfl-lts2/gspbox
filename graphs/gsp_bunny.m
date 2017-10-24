function [ G ] = gsp_bunny()
%GSP_BUNNY Create a graph of the stanford bunny
%   Usage :  G = gsp_bunny();
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_bunny()' creates a graph from the pointcloud of the Stanford Bunny model. 
%
%   Example:::
%
%           G = gsp_bunny();
%           gsp_plot_graph(G);
%
%   References: turk1994zippered

    %Load the point cloud
    P = gsp_pointcloud('bunny');
    
    %Create the graph from the point cloud using an epsilon-neighborhood
    %connectivity
    param.type = 'knn';
    param.rescale = 1;
    param.center = 1;
    
    %Compute it
    G = gsp_nn_graph(double(P), param);
    %Reduce vertex size for plotting
    G.plotting.vertex_size = 10;
end

