function A = createRandRegGraph(vertNum, deg)
% createRegularGraph - creates a simple d-regular undirected graph
% simple = without loops or double edges
% d-reglar = each vertex is adjecent to d edges
%
% input arguments :
%   vertNum - number of vertices
%   deg - the degree of each vertex
%
% output arguments :
%   A - A sparse matrix representation of the graph
%
% algorithm :
% "The pairing model" : create n*d 'half edges'.
% repeat as long as possible: pick a pair of half edges 
%   and if it's legal (doesn't creat a loop nor a double edge)
%   add it to the graph
% reference: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.7957&rep=rep1&type=pdf
%
% Written by Golan Pundak and downloaded from the MATLAB Central File
% Exchange
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/misc/createRandRegGraph.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

n = vertNum;
d = deg;
matIter = 10;

%check parameters
if mod(n*d,2)==1   
    disp('createRandRegGraph input err: n*d must be even!');
    A=[];
    return;
end

%a list of open half-edges
U = repmat(1:n,1,d);

%the graphs adajency matrix
A=sparse(n,n);

edgesTested=0; 
repetition=1;

%continue until a proper graph is formed
while ~isempty(U) && repetition < matIter
    
    edgesTested = edgesTested + 1;

    %print progress
    if mod(edgesTested, 5000)==0 
        fprintf('createRandRegGraph() progress: edges=%d/%d\n', edgesTested, n*d);    
    end

    %chose at random 2 half edges
    i1 = ceil(rand*length(U));
    i2 = ceil(rand*length(U));
    v1 = U(i1);
    v2 = U(i2);

    %check that there are no loops nor parallel edges
    if (v1 == v2) || (A(v1,v2) == 1)
        
        %restart process if needed
        if (edgesTested == n*d)           
            repetition=repetition+1;            
            edgesTested = 0;
            U = repmat(1:n,1,d);
            A = sparse(n,n);
        end
    else
        %add edge to graph
        A(v1, v2)=1;
        A(v2, v1)=1;
        
        %remove used half-edges
        v = sort([i1,i2]);
        U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];
    end
end

%check that A is indeed simple regular graph
msg=isRegularGraph(A);
if ~isempty(msg)    
    disp(msg);
end

%-------------------------------------------------

function msg=isRegularGraph(G)
%is G a simple d-regular graph the function returns []
%otherwise it returns a message describing the problem in G

msg=[];

%check symmetry
if (norm(G-G','fro')>0)
    msg=[msg,' is not symmetric, '];
end

%check parallel edged
if (max(G(:))>1)
    msg=[msg,sprintf(' has %d parallel edges, ',length(find(G(:)>1)) )];
end

%check that d is d-regular
d_vec=sum(G);
if min(d_vec)<d_vec(1) || max(d_vec)>d_vec(1)
    msg=[msg,' not d-regular, '];
end

%check that g doesn't contain any loops
if (norm(diag(G))>0)
    msg=[msg,sprintf(' has %d self loops, ',length(find(diag(G)>0)) )];
end

   



