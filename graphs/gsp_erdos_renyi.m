function [ G ] = gsp_erdos_renyi( N,p,param )
%GSP_ERDOS_RENYI Create a random Erdos Renyi graph
%   Usage:  G = gsp_erdos_renyi( N,p,param );
%           G = gsp_erdos_renyi( N,p );
%
%   Input parameters:
%         N     : Number of nodes
%         p     : Probability of connection of a node with another
%         param : Structure of optional parameter
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_erdos_renyi( N,p,param )' initializes create a graph structure
%   containing the weighted adjacency matrix (G.W), the number of vertices
%   (G.N) for an Erdos Renyi graph. All edge weights are equal to 1. 
% 
%   The Erdos Renyi graph is constructed by connecting nodes randomly. Each 
%   edge is included in the graph with probability p independent from every
%   other edge. 
%
%   *param* a Matlab structure containing the following fields:
%
%     * *param.connected* : flag to force the graph to be connected. By
%       default, it is 0.
%
%     * *param.maxit* : is the maximum number of try to connect the graph. 
%       By default, it is 10.
% 
%     * *param.verbose* : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%
%   Example::
%
%        G = gsp_erdos_renyi(100,0.05)

% Author : Nathanael Perraudin
 

% Optional input arguments
if nargin<3, param=struct; end

if ~isfield(param, 'connected'), param.connected=0; end
if ~isfield(param, 'maxit'), param.maxit=10; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end


% Check if the user didn't put silly parameter
if param.connected
    if N*p<3
        if param.verbose
            fprintf(['The model has very low chance to create a '...
                ,'connected graph. Increase N, p or disable '...
                ,'param.connected \n']);
        end
    end
end

if p>1
   error('The probability p cannot be above 1'); 
end

if p<0
    error('The probability p cannot be negative');
end



% Create the graph structure
G=struct;
G.graph_type='erdos_renyi';
G.connected=param.connected;

% Create the graph
if param.connected
    bool=0;
    iter=1;
    if param.verbose>1
        fprintf('Entering the iteration process\n')
    end
    while ~bool
        
        G.W = sprandsym(N, p)>0;
        G.W(1:N+1:end) = 0;

        bool=gsp_check_connectivity(G.W);
        if iter==param.maxit
            warning('The graph is not strongly connected.');
            break;
        elseif ~bool
            if param.verbose>1
                fprintf('Iteration %i failed. Trying again \n',iter)
            end
            iter=iter+1;
        end
    end
    if param.verbose>1
        if bool
            fprintf('A connected graph has been created in %i iteration\n',iter);
        else
            fprintf('The iterative process has failed.\n');
        end
    end
else
   G.W = sprandsym(N,p)>0;
   G.W(1:N+1:end) = 0;
end
    

% try to create point to be associated with the points of the graph

% [x, y] = draw_dot(G.W);
% 
% 
% G.coords=[x',y'];
% G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];

G.type = 'erdos_renyi';

G = gsp_graph_default_parameters(G);

end
