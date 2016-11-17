% function MOD = compute_modularity(IDX,W)
%
% This computes the modularity MOD of a given partition IDX of a given 
% graph represented by its adjacency matrix W
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

function MOD = compute_modularity(IDX,W)

m = sum(sum(W));
MOD = 0;
COMu = unique(IDX);
for j=1:length(COMu)
    IDXj = find(IDX==COMu(j));
    Ec = sum(sum(W(IDXj,IDXj)));
    Et = sum(sum(W(IDXj,:)));
    if Et>0
        MOD = MOD + Ec/m-(Et/m)^2;
    end
end
