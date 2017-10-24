function [G] = gsp_minnesota(connect)
%GSP_MINNESOTA  Initialize the Minnesota road network
%   Usage:  G = gsp_minnesota();
%           G = gsp_minnesota(connect);
%
%   Input parameters:
%         connect : change the graph so that it is connected (default 1)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_minnesota()' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.limits) of the Minnesota road network from the MatlabBGL library.
%
%   Remark: if connect is set to 1. We adjust the adjacency matrix so that
%   all edge weights are equal to 1, and the graph is connected.  It is the
%   default!
%
%   To get the orinial disconnected graph, use::
%
%           G = gsp_minnesota(connect);
%
%   Example:::
%
%          G = gsp_minnesota();
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   References: gleich

% Author : David I Shuman, Nathanael Perraudin
% Date : 15 March 2014

if nargin < 1
   connect = 1; 
end

Q=load('minnesota.mat');
G.N=size(Q.A,1);
G.coords=Q.xy;
G.plotting.limits=[-98,-89,43,50];

if connect
    % Edit adjacency matrix
    A=Q.A;
    % clean minnesota graph
    A=A-diag(diag(A));
    % missing edge needed to connect graph
    A(349,355)=1;
    A(355,349)=1;
    % change a handful of 2 values back to 1
    A(86,88)=1;
    A(88,86)=1;
    A(345,346)=1;
    A(346,345)=1;
    A(1707,1709)=1;
    A(1709,1707)=1;
    A(2289,2290)=1;
    A(2290,2289)=1;
    G.W=sparse(A);
    G.type = 'minnesota';
else
    G.W=sparse(Q.A);
    G.type = 'minnesota-disconnected';
end

G.plotting.vertex_size=30;

G = gsp_graph_default_parameters(G);

end

