function gsp_plot_edges(G,param)

if nargin<2
    param = struct;
end

if ~isfield(param,'edge_size'), param.edge_size = 0; end


[ki,kj]=find(G.W);
    if G.directed

        if size(G.coords,2)==2
            In  = [G.coords(ki,1),G.coords(ki,2)]; 
            Fin = [G.coords(kj,1),G.coords(kj,2)]; 
            V=Fin-In;
            quiver(G.coords(ki,1),G.coords(ki,2),V(:,1),V(:,2),0,...
                    '-r','LineWidth',G.plotting.edge_width);
        else
            In  = [G.coords(ki,1),G.coords(ki,2),G.coords(ki,3)]; 
            Fin = [G.coords(kj,1),G.coords(kj,2),G.coords(kj,3)]; 
            V=Fin-In;
            quiver3(G.coords(ki,1),G.coords(ki,2),G.coords(ki,3),...
                V(:,1),V(:,2),V(:,3),0,...
                    '-r','LineWidth',G.plotting.edge_width);
        end
    else
        %gplot23D(G.A,G.coords,G.plotting.edge_style,'LineWidth',G.plotting.edge_width,'Color',G.plotting.edge_color); 
        if param.edge_size
            if numel(ki) > 1000
                disp('Plotting edges with different sizes. It may take some time.');
            end
            if size(G.coords,2) == 2
                already_hold = ishold;
                if not(already_hold)
                    hold on
                end
                for ii = 1:numel(ki)
                    plot([G.coords(ki(ii),1)';G.coords(kj(ii),1)'],...
                        [G.coords(ki(ii),2)';G.coords(kj(ii),2)'],...
                        G.plotting.edge_style, 'LineWidth',...
                        G.plotting.edge_width * (3*G.W(ki(ii),kj(ii))+0.01) ,...
                        'Color',G.plotting.edge_color);
%                     G.plotting.edge_width * (2*G.W(ki(ii),kj(ii))+0.01)
                end
                if not(already_hold)
                    hold off;
                end
            else
                    plot3([G.coords(ki,1)';G.coords(kj,1)'],...
                    [G.coords(ki,2)';G.coords(kj,2)'],...
                    [G.coords(ki,3)';G.coords(kj,3)'],...
                    G.plotting.edge_style, 'LineWidth',G.plotting.edge_width,...
                    'Color',G.plotting.edge_color);
            end
        else
            if size(G.coords,2) == 2
                plot([G.coords(ki,1)';G.coords(kj,1)'],...
                    [G.coords(ki,2)';G.coords(kj,2)'],...
                    G.plotting.edge_style, 'LineWidth',G.plotting.edge_width,...
                    'Color',G.plotting.edge_color);
            else
                    plot3([G.coords(ki,1)';G.coords(kj,1)'],...
                    [G.coords(ki,2)';G.coords(kj,2)'],...
                    [G.coords(ki,3)';G.coords(kj,3)'],...
                    G.plotting.edge_style, 'LineWidth',G.plotting.edge_width,...
                    'Color',G.plotting.edge_color);
            end
        end
    end
end