function gsp_plot_graph(G,param)
%GSP_PLOT_GRAPH  Plot a graph in 2D or 3D
%   Usage:  gsp_plot_graph(G);
%           gsp_plot_graph(G,param);
%
%   Input parameters:
%         G          : Graph structure or a cell of graph structures.
%         param      : Optional variable containing additional parameters.
%   Output parameters:
%         none
%
%   'gsp_plot_graph(G)' plots a graph (or multiple graphs) in 2D or 3D, 
%   using the adjacency matrix (G.A), the plotting coordinates (G.coords), 
%   the coordinate limits (G.plotting.limits), the edge width (G.plotting.edge_width), 
%   the edge color (G.plotting.edge_color), the edge style (G.plotting.edge_style), the 
%   vertex size (G.plotting.vertex_size), and the vertex color (G.vertex_color).
%
%   Additional parameters:
%   * *param.show_edges*     : Set to 0 to only draw the vertices. (default G.Ne < 10000 )
%   * *param.num_clusters*   : Number of clusters for a clustered graph.
%   * *param.clusters*       : Cluster identities for a clustered graph.
%   * *param.cluster_colors* : Cluster colors for a clustered graph.
%   * *param.cp*             : Camera position for a 3D graph
%
%   Example:::
%
%          G = gsp_swiss_roll(200);
%          gsp_plot_graph(G);
%
%   See also: gsp_plot_signal gsp_plot_signal_spectral
%


% Author: David I Shuman, Nathanael Perraudin
% Date: 14 March 2014

  
% Handling optional inputs
if nargin<2
   param = struct;
end

if ~numel(G.coords)
    error('There is no coordinates associated to this graph!');
end


if ~isfield(param,'show_edges')  
	param.show_edges=G.Ne<10000;
end

if ~isfield(param,'clear'), param.clear = 1; end
if ~isfield(param,'edge_size'), param.edge_size = 0; end


  
if iscell(G)
    for i=1:length(G)
        figure;
        gsp_plot_graph(G{i},param);
    end
else
    
    if ~isfield(G,'coords')
        error('GSP_PLOT_GRAPH: Cannot plot a graph without coordinate!')
    end
    
    % Clear axes
    if param.clear
        cla;
    end
    G = gsp_graph_default_plotting_parameters(G);
    
    % TODO: To be changed
    if (~isfield(param,'num_clusters') || ...
            ~isfield(param,'clusters') || ...
            ~isfield(param,'cluster_colors') )
        num_clusters=1;
        clusters=ones(G.N,1);
        cluster_colors=G.plotting.vertex_color;
    elseif length(param.cluster_colors)<param.num_clusters 
        error('Not enough cluster colors specified');
    else
        num_clusters=param.num_clusters;
        clusters=param.clusters;
        cluster_colors=param.cluster_colors;
    end

    % Start ploting
    hold on;
    
    
    if param.show_edges
        gsp_plot_edges(G,param);
    end

    for i=1:num_clusters

        if size(G.coords,2)==2
            scatter(G.coords((clusters==i),1),G.coords((clusters==i),2),...
            G.plotting.vertex_size,'MarkerFaceColor',cluster_colors(i, :),...
            'MarkerEdgeColor',G.plotting.vertex_color);

        elseif size(G.coords,2)==3
            scatter3(G.coords((clusters==i),1),G.coords((clusters==i),2),...
                G.coords((clusters==i),3),G.plotting.vertex_size,...
                'MarkerFaceColor',cluster_colors(i, :),'MarkerEdgeColor',...
                G.plotting.vertex_color);
%             xlim([G.plotting.limits(1,1) G.plotting.limits(2,1)]);
%             ylim([G.plotting.limits(1,2) G.plotting.limits(2,2)]);
%             zlim([G.plotting.limits(1,3) G.plotting.limits(2,3)]);
        end
        axis(G.plotting.limits);
    end

    if size(G.coords,2)==3
        if ~isfield(param,'cp')
            cp=[-1.4,-16.9,3.4];
        else 
            cp=param.cp;
        end
        set(gca,'CameraPosition',cp);
    end

    axis off;
    hold off;
end

% TODO: use special presets of the following style for more graphs:
if isfield(G, 'type')
    if strcmp(G.type, 'non_uniform') ||...
       strcmp(G.type, 'sub-non_uniform') ||...
       strcmp(G.type, 'non_uniform_patch') ||...
       strcmp(G.type, 'sub-non_uniform_patch')
        axis equal; axis tight;% axis on;
    end
end

end
