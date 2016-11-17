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
%                     G.plotting.edge_width  (2*G.W(ki(ii),kj(ii))+0.01)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/plotting/gsp_plot_edges.php

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
