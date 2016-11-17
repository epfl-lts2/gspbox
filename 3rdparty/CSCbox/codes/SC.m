% function [IDX_SC, Uk, Dk, time_SC] = SC(G,lap_type)
%
% This is an implementation of the classical Spectral CLustering algorithm 
% using either: 
% - the combinatorial Laplacian if lap_type='combinatorial' (in this case, 
% the features are not normalized);
% - or the normalized Laplacian if lap_type='normalized' (in this case, 
% the features are normalized). 
% See the paper cited below for many references on this algorithm. 
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


function [IDX_SC, Uk, Dk, time_SC] = SC(G,lap_type)

normBOOL=strcmp(lap_type,'normalized');
G = gsp_create_laplacian(G,lap_type);

% eigs
if isfield(G, 'G.Dk'), % if it was already calculated, do not do it again
    Dk=G.Dk;
    Uk=G.Uk;
    time_SC.eigs=G.time_SC_eigs;
else % run eigs
    tic;
    opts.isreal=1;opts.issym=1;opts.maxit=10000;
    [Uk,Dk]=eigs(G.L,G.k,'SA',opts);Dk=diag(Dk);
    time_SC.eigs=toc;
end

if normBOOL % kmeans on Yk
    tic;
    Yk=Uk./repmat(sqrt(sum(Uk.^2,2)),1,G.k);
    IDX_SC = kmeans(Yk, G.k,'Replicates',20);
    time_SC.kmeans=toc;
else % kmeans on Uk
    tic;
    IDX_SC = kmeans(Uk, G.k,'Replicates',20);
    time_SC.kmeans=toc;
end

time_SC.total=time_SC.eigs+time_SC.kmeans;
