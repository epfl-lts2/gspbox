function [G]=gsp_airfoil()
%GSP_AIRFOIL Initialize the airfoil graph
%   Usage:  G=gsp_airfoil();
%
%   Input parameters:
%         none
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_airfoil()' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.limits) of the airfoil mesh graph. All edge weights are equal
%   to 1.
%
%   Example:::
%
%          G = gsp_airfoil();
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   See also: gsp_graph


% Author : David I Shuman, Nathanael Perraudin
% Date : 15 March 2014

load airfoil
%G.N=4253;
A = sparse(i_inds,j_inds,1,4253,4253);
G.W = (A + A')/2;
G.coords=[x,y];
G.plotting.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];

G.type = 'airfoil';

G.plotting.vertex_size=30;

G = gsp_graph_default_parameters(G);


end

