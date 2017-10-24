%FRONTPAGE

close all;
clear;

paramplot.save = 1;
paramplot.savefig = 0;
paramplot.pathfigure = [fileparts(mfilename('fullpath')),'/'];
paramplot.eps = 1;



G = gsp_logo;
G.plotting.vertex_size = 200;
% G.plotting.edge_color = [200/255 136/255.0 204/255.0];
G.plotting.edge_color = [0.5 0.5 0.5];
G.plotting.edge_width = 2;
f = zeros(G.N,1);
f(G.info.idx_g) = -1;
f(G.info.idx_s) = 1;
f(G.info.idx_p) = 0;
f = f + 0.2*randn(G.N,1);

figure;
paramplot.position = [100,100,1200,600];
gsp_plot_signal(G,f)
colorbar off
gsp_plotfig('frontpage',paramplot)

