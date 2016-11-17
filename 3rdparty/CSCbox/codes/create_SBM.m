% function [A,truth]=create_SBM(N,q,c,epsi,com_size)
%
% This creates an SBM graph with
% Inputs:
% - N nodes, 
% - q communites of sizes listed in com_size, 
% - an average degree of c,
% - and a difficulty \epsilon=epsi. 
% Outputs:
% - G the graph structure
%
% Copyright (C) 2016 Nicolas Tremblay, Gilles Puy.
% This file is part of the CSCbox (Compressive Spectral Clustering toolbox)
%
% The CSCbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The CSCbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%     N. Tremblay, G. Puy, R. Gribonval and P. Vandergheynst.
%     Compressive Spectral Clustering.
%     ArXiv e-prints:1602.02018 Feb. 2016.

function [G]=create_SBM(N,q,c,epsi,com_size)


if sum(com_size)~=N
    error('error in create_SBM.m : sum(com_size)~=N!')
end

pin=(q*c)/(N-q+(q-1)*epsi*N);
pout=(q*c*epsi)/(N-q+(q-1)*epsi*N);

rng('shuffle');

% total # of intra community edges is taken from Bernouilli distrib
I=[]; J=[]; truth=zeros(N,1);
for k=1:q
    truth(sum(com_size(1:k-1))+1:sum(com_size(1:k)))=k;
    Numedgesk=my_binornd(com_size(k).*(com_size(k)-1)/2,pin); %normal approx
    while Numedgesk>com_size(k).*(com_size(k)-1)/2
        Numedgesk=my_binornd(com_size(k).*(com_size(k)-1)/2,pin); %normal approx
    end
    edges = randperm(com_size(k).*(com_size(k)-1)/2,Numedgesk);
    %edges=datasample([1:com_size(k).*(com_size(k)-1)./2],Numedgesk,'Replace',false);
    [Inew, Jnew] = ind2sub4up(edges);
    I=[I;Inew+sum(com_size(1:k-1))];
    J=[J;Jnew+sum(com_size(1:k-1))];
end

% for k=1:q-1
%     for kk=k+1:q
%         Numedgesk=my_binornd(com_size(k)*com_size(kk),pout); %normal approx
%         edges = randperm(com_size(k)*com_size(kk),Numedgesk);
%         %Numedgesk=binornd(com_size(k)*com_size(kk),pout);
%         %edges=datasample([1:com_size(k)*com_size(kk)],Numedgesk,'Replace',false);
%         [Inew, Jnew] = ind2sub([com_size(k),com_size(kk)],edges);
%         I=[I;Inew'+sum(com_size(1:k-1))];
%         J=[J;Jnew'+sum(com_size(1:kk-1))];
%     end
% end

for k=1:q-1
    Numedgesk=my_binornd(com_size(k)*sum(com_size(k+1:end)),pout); %normal approx
    edges = randperm(com_size(k)*sum(com_size(k+1:end)),Numedgesk);
    %Numedgesk=binornd(com_size(k)*sum(com_size(k+1:end)),pout);
    %edges=datasample([1:com_size(k)*sum(com_size(k+1:end))],Numedgesk,'Replace',false);
    [Inew, Jnew] = ind2sub([com_size(k),sum(com_size(k+1:end))],edges);
    I=[I;Inew'+sum(com_size(1:k-1))];
    J=[J;Jnew'+sum(com_size(1:k))];
end

A=sparse(I,J,1,N,N);
A=A+A';
G=gsp_graph(A);
G.type='SBM_CSC';
G.truth = truth;
G.k=q;