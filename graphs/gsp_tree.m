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
%   Example:
%
%          G = gsp_tree(4,3);
%          gsp_plot_graph(G,paramplot);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_tree.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781


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


