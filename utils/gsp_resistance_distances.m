function [resistance_distances] = gsp_resistance_distances(G,param)
%GSP_RESISTANCE_DISTANCES : Compute the resitance distances of a graph
%   Usage: rd = gsp_resistance_distances(G);
%          rd = gsp_resistance_distances(L);
%
%   Input parameters:
%       G    : Graph structure or Laplacian matrix (L)
%       param: optional parameters
%   Output parameters:
%       rd   : distance matrix
%
%   This function compute the resistance distances between all vertices in 
%   a graph. The distance between two nodes is defined as the inverse of 
%   the weight matrix. For example the distance matrix::
%
%           dist = [0, 3, 1;...
%                   3, 0, 2;...
%                   1, 2, 0];
%
%   corresponds to the weight matrix::
%
%           W = [0, 1/3, 1/1;...
%                1/3, 0, 1/2;...
%                1/1, 1/2, 0];
%
%   The function will compute the resistance distance following 
%   Kirshoff's law. In the our example it is::
%
%           rd2 = [0, 3/2, 5/6;...
%                  3/2, 0, 4/3;...
%                  5/6, 4/3, 0]
%
%   In matlab, you can reproduce this example using:::
%
%           % The weight 
%           dist = [0, 3, 1;...
%                   3, 0, 2;...
%                   1, 2, 0];
%           % The weight is the inverse of the distance...
%           W = dist.^(-1);
%           % Fix the diagonal
%           W([1,5,9])=0;    
%           G = gsp_graph(W);
%           rd = gsp_resistance_distance(G)
%           % Resitance computed by hand
%           rd2 = [0, 3/2, 5/6;...
%                  3/2, 0, 4/3;...
%                  5/6, 4/3, 0]
%
%   *param* is an optional structure that contains the following field
%
%   * *param.verbose*: display parameter - 0 no log - 1 display warnings
%     (default 1)
%   
%   References: klein1993resistance

% Author: Nathanael Perraudin, David Shuman
% Testing: test_resistance_distance

if nargin < 2
    param = struct;
end

if ~isfield(param,'verbose')
    param.verbose = 1;
end

if isstruct(G)
    % Use the non normalized laplacian
    if ~strcmp(G.lap_type, 'combinatorial')
        G = gsp_create_laplacian(G,'combinatorial');
        if param.verbose
            fprintf(['Compute the combinatorial laplacian ',...
            'for the resitance distance\n']);
        end
    end
    L = G.L;
    if (isfield(G,'U') || isfield (G,'e') )
        if issorted(G.e)
            pseudo=G.U(:,2:G.N)*diag(1./G.e(2:G.N))*G.U(:,2:G.N)';
        else
            [e, si] = sort(G.e,'ascend');
            U = G.U(:,si);
            pseudo=U(:,2:G.N)*diag(1./e(2:G.N))*U(:,2:G.N)';
        end
    else 
        pseudo=pinv(full(L));
    end
else
    L = G;
    pseudo=pinv(full(L));
end

N=size(L,1);
d = diag(pseudo);
resistance_distances = repmat(d,1,N) + repmat(d',N,1) - pseudo - pseudo';


