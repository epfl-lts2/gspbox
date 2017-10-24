function [G]=gsp_low_stretch_tree(k)
%GSP_LOW_STRETCH_TREE  Initialize a low stretch tree
%   Usage:  G = gsp_low_stretch_tree(k);
%           G = gsp_low_stretch_tree();
%
%   Input parameters:
%         k     : 2^k points on each side of the grid of vertices. (default 6)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_create_low_stretch_tree(k)' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), the plotting coordinate limits 
%   (G.limits), and the root of a low stretch tree on a grid of points.
%   There are $2^k$ points on each side of the grid, and therefore $2^{2k}$ 
%   total vertices. The edge weights are all equal to 1.
%
%   Example:::
%
%          G = gsp_low_stretch_tree(3);
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%   See also: gsp_graph


% Author : David I Shuman, Nathanael Perraudin

if nargin < 1
   k = 6; 
end

ii=[2 3 1 1 4 3];
jj=[1 1 2 3 3 4];

XCoords=[1, 2, 1, 2];
YCoords=[1, 1, 2, 2];

for p=2:k
    % create the indicies of the weight matrice
    ii_new=[ii,ii+4^(p-1),ii+2*4^(p-1),ii+3*4^(p-1)];
    ii_new_middle=[4^(p-1),4^(p-1),4^(p-1)+(4^p+2)/3,5/3*4^(p-1)+1/3,4^(p-1)+(4^p+2)/3,3*4^(p-1)+1];
    ii=[ii_new,ii_new_middle];

    jj_new=[jj,jj+4^(p-1),jj+2*4^(p-1),jj+3*4^(p-1)]; 
    jj_new_middle=[5/3*4^(p-1)+1/3,4^(p-1)+(4^p+2)/3,3*4^(p-1)+1,4^(p-1),4^(p-1),4^(p-1)+(4^p+2)/3];
    jj=[jj_new,jj_new_middle];
    
    % Create Coords
    YCoords=repmat(YCoords,1,2);
    YCoords_new=[YCoords,YCoords+2^(p-1)];
    YCoords=YCoords_new;

    XCoords_new=[XCoords,XCoords+2^(p-1)];
    XCoords=repmat(XCoords_new,1,2);
end

G.W=sparse(ii,jj,ones(size(ii)));
G.coords=[XCoords',YCoords'];

G.limits=[0,2^k+1,0,2^k+1];
G.N=(2^k)^2;
G.root=4^(k-1);

G.plotting.edge_width=1.25;
G.plotting.vertex_size=75;

G.type = 'low strech tree';

G = gsp_graph_default_parameters(G);


end
