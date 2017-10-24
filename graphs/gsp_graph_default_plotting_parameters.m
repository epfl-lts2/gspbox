function G = gsp_graph_default_plotting_parameters(G)
%GSP_GRAPH_DEFAULT_PLOTTING_PARAMETERS Default plotting parameters
%   Usage: G = gsp_graph_default_plotting_parameters( G );
%
%   Input parameters
%       G   : Graph 
%   Output parameters
%       G   : Graph
% 
%   This function complete all plotting parameter for a graph.
%
%
%   List of plotting parameters
%   ---------------------------
%
%   
%   *G.plotting*: Plotting parameters
%   * *G.plotting.edge_width*: Width of edges (default 1)
%   * *G.plotting.edge_color*: Color of edges (default [255,88,41]/255 )
%   * *G.plotting.edge_style*: Style of edges (default '-')
%   * *G.plotting.vertex_size*: Size of vertex (default 50)
%   * *G.plotting.vertex_color*: Color of vertex (default 'b')
%   
%          

% Author: Nathanael Perraudin
% Date  : 13.08.2014

if nargin<1
    error('You need a graph to run this function!')
end

if isfield(G,'coords')

    % Assign default plotting parameters
    if ~isfield(G,'plotting')
        G.plotting=struct; 
    end


    if ~isfield(G.plotting,'edge_width')
        G.plotting.edge_width=1; 
    end
    if ~isfield(G.plotting,'edge_color')
    %     G.plotting.edge_color=[255,88,41]/255;
        G.plotting.edge_color=[0.6 0.6 0.6];
    end
    if ~isfield(G.plotting,'edge_style')
        G.plotting.edge_style='-';
    end
    if ~isfield(G.plotting,'vertex_size')
        G.plotting.vertex_size=50;
    end
    if ~isfield(G.plotting,'vertex_color')
        G.plotting.vertex_color='b';
    end
    % if ~isfield(G,'vertex_edge_color')
    %     G.plotting.vertex_edge_color='k';
    % end


    % Limit for plotting
    if ~isfield(G.plotting,'limits'), 
        if isfield(G,'coords')
            if size(G.coords,2)>2
                mrx =  (max(G.coords(:,1)) - min(G.coords(:,1)))/20+eps;
                mry =  (max(G.coords(:,2)) - min(G.coords(:,2)))/20+eps;
                mrz =  (max(G.coords(:,3)) - min(G.coords(:,3)))/20+eps;
                G.plotting.limits = [min(G.coords(:,1))-mrx,...
                    max(G.coords(:,1))+mrx, ...
                    min(G.coords(:,2))-mry, max(G.coords(:,2))+mry, ... 
                    min(G.coords(:,3))-mrz, max(G.coords(:,3))+mrz]; 
            elseif size(G.coords,2)>1
                mrx =  (max(G.coords(:,1)) - min(G.coords(:,1)))/20+eps;
                mry =  (max(G.coords(:,2)) - min(G.coords(:,2)))/20+eps;
                G.plotting.limits = [min(G.coords(:,1))-mrx,...
                    max(G.coords(:,1))+mrx, min(G.coords(:,2))-mry,...
                    max(G.coords(:,2))+mry]; 
            else
                mrx =  (max(G.coords) - min(G.coords))/20;
                G.plotting.limits = [min(G.coords)-mrx, max(G.coords)+mrx]; 
            end

        else
            G.plotting.limits = [0, 1, 0,1,0, 1]; 
        end
    end
else
    G.plotting = struct;
end

end
