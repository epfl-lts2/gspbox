function gsp_plot_jtv_signal(G,signal,param)
%GSP_PLOT_JTV_SIGNAL  Plot a time-vertex signal in 2D or 3D
%   Usage:  gsp_plot_jtv_signal(G,signal);
%           gsp_plot_jtv_signal(G,signal,param);
%
%   Input parameters:
%         G          : Time-Vertex graph structure.
%         signal     : Time-Vertex signal.
%         param      : Optional variable containing additional parameters.
%   Output parameters:
%         none
%
%   'gsp_plot_jtv_signal(G,signal)' plots a time-vertex signal in 2D or 3D,
%   using the adjacency matrix (G.A), the plotting coordinates (G.coords),
%   the coordinate limits (G.plotting.limits), the edge width
%   (G.plotting.edge_width), the edge color (G.plotting.edge_color), the
%   edge style (G.plotting.edge_style), and the vertex size
%   (G.vertex_size). If the figure loses the focus the visualization will
%   stop
%
%   Example:::
%
%          G = gsp_sensor();
%          G = gsp_jtv_graph(G,50);
%          f = sin(2*pi/G.N*(1:G.N)'*[1:G.jtv.T]);
%          gsp_plot_jtv_signal(G,f)
%
%   Additional parameters
%   ---------------------
%
%   * *param.step*       : Time step (default round(log(param.T)))
%   * *param.T*          : Final time (use if you want to visualize only the first param.T samples)
%   * *param.show_edges* : Set to 0 to only draw the vertices. (default 0, really slow )
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
%   See also: gsp_plot_signal
%


if nargin < 3
    param = struct;
end


global brk;brk=1;
hFig=figure;
% Get the underlying Java reference
warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
jFig = get(hFig, 'JavaFrame');

jAxis = jFig.getAxisComponent;


if ~isfield(param,'show_edges'), param.show_edges = 0; end
if ~isfield(param,'T'), param.T = size(signal,2); end
if ~isfield(param,'step'), param.step = round(log(param.T)); end
if ~isfield(param,'colorbar'), param.colorbar = 0; end
if ~isfield(param,'climits'), param.climits = [min(signal(:)) max(signal(:))]/2; end
if ~isfield(param,'speed'), param.speed = 0.1; end

step = param.step;
speed = param.speed;
T = param.T;

switch G.type
    case {'ring','path'}
        f=signal(:,1);
        plot(f)
        ylim([min(signal(:)), max(signal(:))])
        title(1)
        drawnow
        set(jAxis.getComponent(1),'FocusLostCallback',{@Lostfocus,hFig});
        ii=2;
        while (ii<T) && brk
            f=signal(:,ii);
            plot(f)
            ylim([min(signal(:)), max(signal(:))])
            title(ii)
            drawnow
            pause(speed)
            ii = ii + step;
        end
    case '2d-grid'
        f=reshape(signal(:,1),sqrt(G.N),[]);
        imagesc(f)
        caxis([min(signal(:)), max(signal(:))])
        if param.colorbar
            colorbar
        end
        title(1)
        drawnow
        set(jAxis.getComponent(1),'FocusLostCallback',{@Lostfocus,hFig});
        ii=2;
        while (ii<T) && brk
            f=reshape(signal(:,ii),sqrt(G.N),[]);
            imagesc(f)
            caxis([min(signal(:)), max(signal(:))])
            if param.colorbar
                colorbar
            end
            title(ii)
            drawnow
            pause(speed)
            ii = ii + step;
        end
    otherwise
        
        f=signal(:,1);
        figure(hFig)
        gsp_plot_signal(G,f,param)
        title(1)
        drawnow
        set(jAxis.getComponent(1),'FocusLostCallback',{@Lostfocus,hFig});
        ii=2;
        while ii<T && brk
            f=signal(:,ii);
            
            gsp_plot_signal(G,f,param)
            title(ii)
            drawnow
            pause(speed)
            ii = ii + step;
        end
end

set(jAxis.getComponent(1),'FocusLostCallback',{});



    function Lostfocus(jAxis, ~, hFig)
        
        disp('Visualization Paused: CTRL + C to interrupt or any key to continue...')
        set(jAxis.getComponent,'FocusLostCallback',{});
        brk=0;
        pause
        set(jAxis.getComponent,'FocusLostCallback',{@Lostfocus,hFig});
        brk=1;
    end


end


