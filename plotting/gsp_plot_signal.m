function gsp_plot_signal(G,signal,param)
%GSP_PLOT_SIGNAL  Plot a graph signal in 2D or 3D
%   Usage:  gsp_plot_signal(G,signal);
%           gsp_plot_signal(G,signal,param);
%
%   Input parameters:
%         G          : Graph structure.
%         signal     : Graph signal.
%         param      : Optional variable containing additional parameters.
%   Output parameters:
%         none
%
%   'gsp_plot_signal(G,f)' plots a graph signal in 2D or 3D, using the adjacency 
%   matrix (G.A), the plotting coordinates (G.coords), the 
%   coordinate limits (G.plotting.limits), the edge width (G.plotting.edge_width), the 
%   edge color (G.plotting.edge_color), the edge style (G.plotting.edge_style), and the 
%   vertex size (G.vertex_size).
%
%   Example:::
%
%          G = gsp_ring(15);
%          f = sin((1:15)*2*pi/15);
%          gsp_plot_signal(G,f)
%
%   Additional parameters
%   ---------------------
%   * *param.show_edges* : Set to 0 to only draw the vertices. (default G.Ne < 10000 )
%   * *param.cp*         : Camera position for a 3D graph.
%   * *param.vertex_size*: Size of circle representing each signal component.
%   * *param.colorbar*   : Set to 0 to not show the colorbar
%   * *param.climits*    : Limits of the colorbar
%   * *param.vertex_highlight*: Vector of indices of vertices to be highlighted
%   * *param.bar*        : 1 Display bar for the graph. 0 Display color
%                          points. (default 0);
%   * *param.bar_width*  : Width of the bar (default 1)
%   * *param.clear*      : Clear axes (default 1)
%
%   See also: gsp_plot_graph gsp_plot_signal_spectral
%


% Author :  Nathanael Perraudin, David I Shuman
% Testing: test_plotting

  

% Handling optional inputs
if nargin < 3
   param = struct;
end


if ~numel(G.coords)
    error('There is no coordinates associated to this graph!');
end

if sum(abs(imag(signal(:))))>1e-10
   error('GSP_PLOT_SIGNAL: can not display complex signal') 
end

signal = real(signal);

if ~isfield(param,'show_edges'), param.show_edges = G.Ne < 5000; end
if ~isfield(param,'bar'), param.bar = 0; end
if ~isfield(param,'clear'), param.clear = 1; end
if ~isfield(param,'bar_width'), param.bar_width = 1; end
if ~isfield(param,'vertex_highlight'), param.vertex_highlight = 0; end
if ~isfield(param,'vertex_size'), 
    if ~isfield(G.plotting,'vertex_size')
        param.vertex_size=500; 
    else
        param.vertex_size=G.plotting.vertex_size*10;
    end
else
    param.vertex_size = param.vertex_size*10; 
end

if ~isfield(param,'colorbar'), param.colorbar = 1; end
if ~isfield(param,'cp')
    if isfield(G.plotting, 'cp')
       param.cp = G.plotting.cp; 
    else
       param.cp=[-6   -3  160]; 
    end
    
end

if ~isfield(G,'coords')
    error('GSP_PLOT_SIGNAL: Cannot plot a graph without coordinate!')
end

% Clear axes
if param.clear
    cla;
end
G = gsp_graph_default_plotting_parameters(G);

hold on;


if param.show_edges
    gsp_plot_edges(G,param);
end

if size(G.coords,2) == 2
    if param.bar
        ind = find(signal < 0);
        plot3([G.coords(ind,1)'; G.coords(ind,1)'],...
            [G.coords(ind,2)'; G.coords(ind,2)'],...
            [zeros(1,length(ind)); reshape(signal(ind),1,[])],...
            G.plotting.edge_style, 'LineWidth',param.bar_width,'color','k');
        ind = find(signal >= 0);
        plot3([G.coords(ind,1)'; G.coords(ind,1)'],...
            [G.coords(ind,2)'; G.coords(ind,2)'],...
            [zeros(1,length(ind)); reshape(signal(ind),1,[])],...
            G.plotting.edge_style, 'LineWidth',param.bar_width,'color','b');
        if any(param.vertex_highlight > 0)
            vh = param.vertex_highlight;
            plot3([G.coords(vh,1)'; G.coords(vh,1)'],...
                [G.coords(vh,2)'; G.coords(vh,2)'],...
                [zeros(1,length(vh)); reshape(signal(vh),1,[])],...
                G.plotting.edge_style, 'LineWidth',2*param.bar_width,'color','m');

        end   
    else
        scatter(G.coords(:,1),G.coords(:,2), ...
            param.vertex_size, signal, '.');
        if any(param.vertex_highlight > 0)
            vh = param.vertex_highlight;
            scatter(G.coords(vh,1),G.coords(vh,2), ...
                param.vertex_size/3, 'ok','LineWidth',3);

        end   
    end


else %if size(G.coords,2) == 3
    scatter3(G.coords(:,1),G.coords(:,2),G.coords(:,3),...
                    param.vertex_size,signal,'.');
    if any(param.vertex_highlight > 0)
        vh = param.vertex_highlight;
        scatter3(G.coords(vh,1),G.coords(vh,2),G.coords(vh,3), ...
            param.vertex_size/3, 'ok');

    end             

end


axis(G.plotting.limits);

if size(G.coords,2)==3 || param.bar
    set(gca,'CameraPosition',param.cp);
end

if ~param.bar
    if isfield(param, 'climits')
        caxis(param.climits);
    %else
        % no need to do anything, the following happens anyway!
        %set(gca,'CLimMode','auto');
    end

    if param.colorbar
        colorbar;
    end
end
colormap(jet)
axis off;
hold off;

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
