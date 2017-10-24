function [G]=gsp_torus(N,M)
%GSP_TORUS  Initialize a 2 dimentional grid graph
%   Usage:  G=gsp_path(N);
%
%   Input parameters:
%         N     : Number of vertices along the first dimention (default 16)
%         M     : Number of vertices along the second dimention (default N)
%   Output parameters:
%         G     : Graph structure.
%
%   The 2dring graph correspond the  graph used for the DFT2. 
%
%   Example:::
%
%          G = gsp_torus(16,20);
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


% Create weighted adjancency matrix
K = 2*N;
J = 2*M;
i_inds = zeros(K*M+J*N,1);
j_inds = zeros(K*M+J*N,1);
for ii = 1:M
    i_inds((ii-1)*K+(1:K)) = (ii-1)*N+[N,1:(N-1),1:N]';
    j_inds((ii-1)*K+(1:K)) = (ii-1)*N+[1:N,N,1:(N-1)]';
end

for ii = 1:M-1
    i_inds(K*M+(ii-1)*2*N+(1:2*N)) = [((ii-1)*N+(1:N)),(ii*N+(1:N))]';
    j_inds(K*M+(ii-1)*2*N+(1:2*N)) = [(ii*N+(1:N)),((ii-1)*N+(1:N))]'; 
end
i_inds(K*M+(M-1)*2*N+(1:2*N)) = [(1:N),((M-1)*N+(1:N))]';
j_inds(K*M+(M-1)*2*N+(1:2*N)) = [((M-1)*N+(1:N)),(1:N)]';


G.W = sparse(i_inds,j_inds,ones(K*M+J*N,1),N*M,N*M);

% Create coordinates
T = 1.5 + sin((0:M-1)*2*pi/M);
U = cos((0:M-1)*(2*pi)/M);
G.coords = [reshape((cos((0:N-1)*(2*pi)/N)'*T),M*N,1),...
        reshape(sin((0:N-1)*(2*pi)/N)'*T,M*N,1),...
        reshape(repmat(U,N,1),M*N,1)];
G.plotting.limits = [-2.5, 2.5, -2.5, 2.5, -2.5, 2.5];

G.type = 'torus';
G.plotting.vertex_size = 30;
G.directed = 0;

G = gsp_graph_default_parameters(G);

end
