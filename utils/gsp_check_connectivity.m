function [ bool, in,out ] = gsp_check_connectivity( G,param )
%GSP_CHECK_CONNECTIVITY Check if the graph G is aperiodic strongly connected
%   Usage: bool=gsp_check_connectivity( G );
%          bool=gsp_check_connectivity( L );
%          bool=gsp_check_connectivity( W );
%          [bool,in,out]=gsp_check_connectivity( ... );
%
%   Input parameters:
%       G,W,L: Graph, Laplacian matrix or Weight martrix
%       param: Optional parameters
%
%   Output parameters
%       bool: Boolean
%       in  : Nodes without any in connections
%       out : Nodes without any out connections
%
%   Test if each node have at least one in connection and one out
%   connection. If this simple test give good results, the function compute
%   the perron vector of G and test it. It might take some time.
%
%   *param* is an optional structure that contains the following field
%
%   * *param.verbose*: display parameter - 0 no log - 1 display the errors
%   

% Date: 6 june 2013
% Author: Nathanael Perraudin

%TODO: Use a clever method

% Handle Input parameters
if nargin<1, error('Not enought inputs parameters!'); end

if nargin<2, param=struct; end
    
if ~isfield(param, 'verbose'), param.verbose = 0; end


if isstruct(G)
    % If the graph is undirected, use the other function...
    if G.directed == 0
       [bool, in] = gsp_check_connectivity_undirected(G,param);
       out = in;
       return;
    end
    A=G.W;
else
    A=G;
end

if ~gsp_isdirected(A)
     [bool, in] = gsp_check_connectivity_undirected(G,param);
      out = in;
     return;
end


warning('This code is really bad!')

% Number of vertex
N=length(A);

% Remove the diagonal
A=A-diag(diag(A));

% Check the connecivity -- simple
    bool=~boolean((sum(1.-(sum(A,1)>0))+sum(1.-(sum(A,2)>0))));
    in=find(1.-(sum(A,1)>0));
    out=(find(1.-(sum(A,2)>0)))';
if param.verbose
    fprintf('   ---   Test if the graph is strongly connected   ---\n');
end
    
if bool
% Check the connectivity -- harder

    % Compute the Probablility matrix
    P=A./repmat(sum(A,2),1,N);

    % Compute the perron vector of P
    [phi,max_eig_P] = eigs(P',2);
    % test if max_eig_P==1
    if abs(max_eig_P(1,1)-1)>10e3*eps;
        bool=0;
        if param.verbose
                fprintf('    The maximum eigenvalue of P is not 1. \n');
        end
    else
        if param.verbose
                fprintf('    The maximum eigenvalue of P is 1. \n');
        end
    end
    phi=phi(:,1);
    % Test if the perron vector is positive
    if sum( phi)<0; 
        phi=-phi;
    end

    if sum(phi<=10e3*eps)
        bool=0;
        if param.verbose
                fprintf('    Null of negative entry in the perron vector. \n');
        end
    else
        if param.verbose
                fprintf('    Stricly positive perron vector. \n');
        end
    end
    
    % see code from undirected for comments about this line
    max_eig_P = eigs(speye(size(P))+P,2);
    if abs(max_eig_P(2,2)-2)<10e3*eps;
        bool=0;
        if param.verbose
                fprintf('    Second eigenvalue of P is 1. \n');
        end
    else
        if param.verbose
                fprintf('    Only one unit eigenvalue of P.\n');
        end
    end
    

    
elseif param.verbose;
    fprintf('     Not every node has an input connection and an output connection! \n');
    fprintf('         No in connections for nodes: %i\n',in );
    fprintf('         No out connections for nodes: %i\n', out);
end

if param.verbose
    fprintf('   ---   End of tests   ---\n');   
    if bool
        fprintf('    The graph is aperiodic strongly connected.\n')
    else
        fprintf('    The graph is not aperiodic stronly connected.\n');
    end
end

end

