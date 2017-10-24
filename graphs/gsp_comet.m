function [G] = gsp_comet(N,k)
%GSP_COMET Initialize a comet graph
%   Usage:  G = gsp_comet(N,k);
%
%   Input parameters:
%         N     : Number of vertices. (default 32)
%         k     : Degree of center vertex. (default 12)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_comet(N,k)' initializes the commet graph. The commet graph is a
%   simple path graph with a star od degree $k$ at its end.
%
%   Example:::
%
%          G = gsp_comet(16,8);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_graph, gsp_ring, gsp_path

% Author : David I Shuman, Nathanael Perraudin
% Date : 15 March 2015

if nargin < 2
    k = 12;
end

if nargin < 1
    N = 32;
end


% Create weighted adjancency matrix
i_inds=[ones(1,k),2:k+1,k+1:N-1,k+2:N];
j_inds=[2:k+1,ones(1,k),k+2:N,k+1:N-1];
G.W=sparse(i_inds,j_inds,ones(1,length(i_inds)),N,N);

% Create coordinates
G.coords=zeros(N,2);
inds=(2:k+1)';
G.coords(2:k+1,:)=[cos((inds-1)*2*pi/k),sin((inds-1)*2*pi/k)];
G.coords(k+2:end,1)=(2:N-k)';
G.plotting.limits=[-2,max(G.coords(:,1)),min(G.coords(:,2)),max(G.coords(:,2))];

G.type = 'comet';

G = gsp_graph_default_parameters(G);

end
