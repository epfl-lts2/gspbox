function LG = gsp_line_graph(G,param)
%GSP_LINE_GRAPH Create Line Graph (or edge-to-vertex dual graph) of graph G
%   Usage: LG = gsp_line_graph(G);
%          LG = gsp_line_graph(G,param);
%
%   Input parameters:
%       G        : Graph structure
%       param    : Structure of optional parameter
%   Output parameters:
%       LG       : Graph structure
%
%  
%   Use *param.coords*  to compute coordinates using halfway point of
%   coordinates of graph G (if available)
 
% Author: Francesco Grassi
% Date  : July 2016
 
 
if nargin<2
    param=struct;
end
 
if ~isfield(G,'B'); G = gsp_incidence(G); end;
if ~isfield(param,'coords'); param.coords=isfield(G,'coords'); end;
 
 
%Line Graph adjacency matrix
W_LG = G.B*G.B';
W_LG = W_LG - diag(diag(W_LG));
 
LG = gsp_graph(W_LG);
 
 
if param.coords
    % coordinates
    halfway = @(x) squareform(bsxfun(@plus,x,x')-diag(2*x))/2;
    ndim = size(G.coords,2);
    coords = zeros(G.Ne,ndim);
    
    %small loop over the graph dimension (2 or 3 iterations)
    for n=1:ndim
        diff_coord=halfway(G.coords(:,n));
        diff_coord = diff_coord(squareform(G.A)>0);
        coords(:,n) = diff_coord;
    end
    
    LG.coords = coords;
    
end
 
 
 
 
 


