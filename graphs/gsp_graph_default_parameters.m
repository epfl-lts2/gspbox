function G = gsp_graph_default_parameters(G)
%GSP_GRAPH_DEFAULT_PARAMETERS load default parameters for graphs
%   Usage: G = gsp_graph_default_parameters(G);
%          G = gsp_graph_default_parameters();
%
%   Input parameters
%       G   : Graph (Optional)
%   Output parameters
%       G   : Graph
% 
%   This function will fill a graph with all missing parameters such that
%   it is compabatible with all functions of the GSPBox. If you create a
%   graph manually, you need to set only the weight matrix *W*. If you have
%   some coordonate, you can also set *G.coords*. *G.coords* is a $N$ x $2$
%   or a $N$ x $3$ matrix with each columns beeing the coordonates in each
%   dimention. Finally, we recommend to set the fiel *G.type* with a name
%   that suits your graph.
%
%   Example::
%          
%          W = rand(30);
%          W = (W + W')/2;
%          G.W = W - diag(diag(W));
%          G = gsp_graph_default_parameters(G)
%
%   This function can be used to update the weights of your graph. It will
%   recompute the Laplacian operator. Warning this function does not
%   perform any change to the Fourier basis::
%
%          G.W = Wnew;
%          G = gsp_graph_default_parameters(G);
%
%
%   List of parameters of the graph structure
%   -----------------------------------------
%
%   By default, a graph structure in the GSPbox contains the following
%   parameters:
%   
%   * *G.W*: Weight matrix (empty by default)
%   * *G.A*: Adacency matrix (constructed with *W*)
%   * *G.N*: Number of nodes (`size(W,1)`)
%   * *G.type*: Type of graph ('unknown' by default)
%   * *G.directed*: 1 if the graph is directed, 0 if not
%   * *G.lap_type*: Laplacian type (default 'combinatorial') See the
%     function |gsp_create_laplacian| for a exhaustive list of the
%     available laplacians.
%   * *G.d*: Degree vector (Computed with *G.W*)
%   * *G.Ne*: Number of edges
%   * *G.coords*: Coordinates of the vertices (default (0,0) )
%   * *G.plotting*: Plotting parameters
%     * *G.plotting.edge_width*: Width of edges (default 1)
%     * *G.plotting.edge_color*: Color of edges (default [255,88,41]/255 )
%     * *G.plotting.edge_style*: Style of edges (default '-')
%     * *G.plotting.vertex_size*: Size of vertex (default 50)
%     * *G.plotting.vertex_color*: Color of vertex (default 'b')
%   
%   Remark: There is redudancy between $A$, $W$, $L$. However, the GSPBox
%   is done in matlab and is not suppose yet to scale to graph sufficiently
%   large that your matlab have memory problem. However, this will be most
%   likely change for milestone 1.0.0. If you do have a urgent need to
%   overcome this problem, please contact the devolper team.
%          

% Author: Nathanael Perraudin
% Date  : 09.12.2013

if nargin<1
    G=struct;
end

if ~isfield(G,'W') % Weight matrix
    G.W = sparse(0);

end

G.A=sparse(G.W>0);
G.N = size(G.W,1);



% Type of graph
if ~isfield(G,'type')
    G.type='unknown'; 
end; 

% if ~isfield(G,'directed')
    G.directed = gsp_isdirected(G); 
% end; 

if ~isfield(G,'hypergraph')
    G.hypergraph = 0; 
end; 

% Create the graph Laplacian
if ~isfield(G,'lap_type')
    G.lap_type='combinatorial';
end

% if ~isfield(G,'L')
    G = gsp_create_laplacian(G);
% end

G.d = full(sum(G.W,2));

% Number of edges
if G.directed
    G.Ne = nnz(G.W);
else
    G.Ne = nnz(G.W)/2;
end

if ~isfield(G,'coords') % Coordonates
    G.coords = []; 
end

G = gsp_graph_default_plotting_parameters(G);

end