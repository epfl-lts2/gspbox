function [G]=gsp_modified_path(W)
%GSP_MODIFIED_PATH  Initialize a modified path graph
%   Usage:  G = gsp_modified_path(W);
%
%   Input parameters:
%         W     : Vector of weights.
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_modified_path(N)' initializes a path graph structure. The node are
%   connected like a path but with the weight given in W. The number of
%   node is the length of W + 1.
%
%   Example:::
%
%          W = ones(15,1);
%          W(8) = 0.1;
%          G = gsp_modified_path(W);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_path
%

% Author: Nathanael Perraudin
% Date: 15 March 2014

N=length(W)+1;
G.N = N;
% Create weighted adjancency matrix
i_inds = [1:N-1,2:N];
j_inds = [2:N,1:N-1];
G.W = sparse(i_inds,j_inds,[W(:);W(:)],N,N);

% Create coordinates
d = zeros(N,1);
for ii=2:N
   d(ii) = d(ii-1)+log(1+1./W(ii-1))/log(2); 
end

G.coords = [d(:), zeros(N,1)];
G.limits = [-max(d)*0.05, max(d)*1.05, -1, 1];

G.type = 'modified path';

G = gsp_graph_default_parameters(G);

end
