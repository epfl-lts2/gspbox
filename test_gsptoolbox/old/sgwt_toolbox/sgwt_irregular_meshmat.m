% sgwt_irregular_meshmat : Adjacency matrix from irregular domain mask 
%
% function A = sgwt_irregular_meshmat(mask)
%
% Computes the adjaceny matrix of graph for given 2-d irregular
% domain. Vertices of graph correspond to nonzero elements of
% mask. Edges in graph connect to (up to) 4 nearest neighbors.
% 
% Inputs : 
%
% mask - binary map of desired 2-d irregular domain
%
% Outputs:
% A - adjacency matrix of graph for irregular domain

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function A =sgwt_irregular_meshmat(mask)
ind=nan(size(mask));
ind(logical(mask))=1:nnz(mask);
N=nnz(mask);

% there will be, at most, 2*N edges
% so 4*N nonzero elements in A

% generate list of edges
% i j 1
% whenever vertex i connects to vertex j
i=zeros(4*N,1);
j=zeros(4*N,1);

% Create array of indices
% ni{k} are pixels that have neighbor type k
% nj{k} are the inidices of pixels of the corresponding neighbor
%
% k=1 'top' k=2 'right' k=3 'bottom' k=4 'left'
offset_list={[1 0],[0 1],[-1 0],[0 -1]};
for k=1:numel(offset_list)
    offset=offset_list{k};
    nmask=shift(mask,offset);
    nind=shift(ind,offset);
    hnm =mask & shift(mask,offset);
    % hnm "has neighbor mask" is one for pixels that have the neighbor
    % defined by offset, zero for pixels that do not have such a neighbor
    ni{k}=ind(hnm);
    nj{k}=nind(hnm);
end
i=[ni{1};ni{2};ni{3};ni{4}]; 
j=[nj{1};nj{2};nj{3};nj{4}];

A=sparse(i,j,ones(size(i)));
