function [G]=gsp_2dgrid(N,M)
%GSP_2dgrid  Initialize a 2 dimentional grid graph
%   Usage:  G=gsp_path(N);
%
%   Input parameters:
%         N     : Number of vertices along the first dimention (default 16)
%         M     : Number of vertices along the second dimention (default N)
%   Output parameters:
%         G     : Graph structure.
%
%   The 2d grid graph correspond the  graph used for the DCT2????
%
%   Example:::
%
%          G = gsp_2dgrid(16);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_path
%
%   References: strang1999discrete

% Author: Nathanael Perraudin
% Date: 15 March 2014
% Testing: test_graph

if nargin <1
   N = 16; 
end

if nargin < 2
   M = N; 
end

G.N=N*M;

% Create weighted adjancency matrix
K = 2*(N-1);
J = 2*(M-1);
i_inds = zeros(K*M+J*N,1);
j_inds = zeros(K*M+J*N,1);
for ii = 1:M
    i_inds((ii-1)*K+(1:K)) = (ii-1)*N+[1:N-1,2:N]';
    j_inds((ii-1)*K+(1:K)) = (ii-1)*N+[2:N,1:N-1]';
end

for ii = 1:M-1
    i_inds(K*M+(ii-1)*2*N+(1:2*N)) = [((ii-1)*N+(1:N)),(ii*N+(1:N))]';
    j_inds(K*M+(ii-1)*2*N+(1:2*N)) = [(ii*N+(1:N)),((ii-1)*N+(1:N))]'; 
end


G.W = sparse(i_inds,j_inds,ones(K*M+J*N,1),N*M,N*M);

% Create coordinates
G.coords = [repmat((0:(N-1))'/N,M,1),...
        reshape(repmat((0:(M-1))/M,N,1),M*N,1) ];
G.plotting.limits = [-1/N, 1+1/N, -1/M, 1+1/M];

G.type = '2d-grid';
G.plotting.vertex_size = 30;

G = gsp_graph_default_parameters(G);

end
