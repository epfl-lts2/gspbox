function [ G,N ] = gsp_tree( h,p )
%GSP_SENSOR Create a random sensor graph
%   Usage: G = gsp_sensor( N );
%          G = gsp_sensor( );
%          G = gsp_sensor( N,param );
%
%   Input parameters
%       h       : Height (default 4)
%       p       : Number of children (default 3)
%   Output parameters
%       G       : Graph structure
%
%   This function creates a tree graph.
%
%   Example:::
%
%          G = gsp_tree(4,3);
%          gsp_plot_graph(G,paramplot);
%


% Date: 12 september 2015

if nargin <2
    p = 3;
end

if nargin <1
    h = 4;
end

%Weights
N=sum(p.^(0:h));%number of nodes
W=zeros(N,N);
%to be changed into sparse matrix
%W is symmetric with ones on the diagonal
Max=sum(p.^[0:h-1]);

for i=1:N
    if i<=Max
        idx1=i+(p-1)*(i-1)+1;
        idx2=idx1+(p-1);
        W(i,idx1:idx2)=1;
        W(idx1:idx2,i)=1;%symmetric
        W(i,i)=1;
    end
end
G.W=sparse(W);

%Coordinates
G.coords=zeros(N,2);
nodes(1)=0;%root
inc=1;
for ii=2:p:N
    nodes(ii)=inc;
        for jj=1:p-1
            nodes(ii+jj)=inc;
        end
    inc=inc+1;
end
[x,y]=treelayout(nodes);
G.coords(:,1)=x;
G.coords(:,2)=y;

%Plotting limits
% num=1;
% G.plotting.limits=[-num,num,-num,num];

%Type
G.type='Tree';

%Default graph parameters
G = gsp_graph_default_parameters(G);
end

