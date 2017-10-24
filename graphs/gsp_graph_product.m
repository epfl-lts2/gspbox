function G = gsp_graph_product(G1,G2,param)
%GSP_GRAPH_PRODUCT Compute graph product between G1 and G2
%   Usage: G = gsp_graph_product(G1,G2, param)
%
%   Input parameters:
%       G1       : Graph 1
%       G2       : Graph 2
%       param    : Struct of parameters
%   Output parameters:
%       G        : Graph product
%
%   This function computes the graph product between two graphs according to rule specified in the parameters structure.
%   Factors graph structs are added to the struct of G.
%  
%   Additional parameters
%   ---------------------
%
%   * *param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%   * *param.rule*   : Graph product rule:
%     + cartesian (default)
%     + kronecker
%     + direct_sum
%     + strong
%

% Author:  Francesco Grassi
% Date: July 2016

if nargin<3
    param = struct;
end

if ~isfield(param,'rule'), param.rule='cartesian';end

if ~isstruct(G1)
    W1 = G1;
else
    W1 = G1.W;   
end

if ~isstruct(G2)
    W2 = G2;
else
    W2 = G2.W;
end



switch param.rule
        
    case 'direct_sum'
        W = blkdiag(W1,W2);
        G = gsp_graph(W);
        
        C1=normc(G1.coords);
        C2=normc(G2.coords);
        C2(:,1)=C2(:,1)+max(abs(C1(:,1)))+.1;
        
        G.coords = [C1;C2];
        
        G.type = 'direct_sum';
        
    
    case 'cartesian'
        W = gsp_cartesian(W1,W2);
        G = gsp_graph(W);
        
        G.type = 'cartesian';

    case {'kronecker','direct'}
        W = kron(W1,W2);
        G = gsp_graph(W);
        G.type = 'kronecker';

    case 'strong'
        W = gsp_strong(W1,W2);
        G = gsp_graph(W);
        G.type = 'strong';   
        
    otherwise
        error('Unknown product rule. You can choose between Direct_Sum, Cartesian, Kronecker, Strong');
end


G.Gf = {G1,G2};


end

