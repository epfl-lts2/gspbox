function [G]=gsp_ring(N,k)
%GSP_RING  Initialize a ring graph
%   Usage:  G = gsp_ring(N);
%           G = gsp_ring(N,k);
%           G = gsp_ring();
%
%   Input parameters:
%         N     : Number of vertices. (default 64)
%         k     : Number of neighbors in each direction (default 1)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_ring(N)' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.coord_limits) of a ring graph with N vertices. Each vertex in the 
%   ring has $2k$ neighbors (maximum value of $k$ is $N/2$). The edge 
%   weights are all equal to 1.
%
%   Example:::
%
%          G = gsp_ring(64);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%

% Author : David I Shuman, Nathanael Perraudin

if nargin < 1
   N = 64; 
end

if nargin < 2
    k = 1;
end

G.N=N;

if k>N/2
    error('Too many neighbors requested');
end

% Create weighted adjancency matrix
if k==N/2
    num_edges=N*(k-1)+N/2;
else
    num_edges=N*k;
end
i_inds=zeros(1,2*num_edges);
j_inds=zeros(1,2*num_edges);

all_inds=1:N;
for i=1:min(k,floor((N-1)/2))
   i_inds((i-1)*2*N+1:(i-1)*2*N+N)=all_inds;
   j_inds((i-1)*2*N+1:(i-1)*2*N+N)=1+mod(all_inds-1+i,N);
   i_inds((i-1)*2*N+N+1:i*2*N)=1+mod(all_inds-1+i,N);
   j_inds((i-1)*2*N+N+1:i*2*N)=all_inds;
end

if k==N/2
   i_inds(2*N*(k-1)+1:2*N*(k-1)+N)=all_inds;
   j_inds(2*N*(k-1)+1:2*N*(k-1)+N)=1+mod(all_inds-1+k,N);
end

G.W=sparse(i_inds,j_inds,ones(1,length(i_inds)),N,N);

%TODO: rewrite G.W without for loops

% Create coordinates
G.coords=[(cos((0:N-1)*(2*pi)/N))',(sin((0:N-1)*(2*pi)/N))'];
G.plotting.limits=[-1,1,-1,1];

if k==1 
    G.type = 'ring';
else
    G.type = 'k-ring';
end

G = gsp_graph_default_parameters(G);

end
