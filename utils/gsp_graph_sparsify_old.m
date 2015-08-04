function [ sparserL , sparserA, sparserW ] = gsp_graph_sparsify_old(L,epsilon)
% Spielman-Srivastava - epsilon should be between 1/sqrt(N) and 1
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_graph_sparsify_old.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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
N=size(L,1);

if ( (epsilon <= 1/sqrt(N)) || (epsilon >1) )
    error('Epsilon out of required range');
end

resistance_distances = gsp_compute_resistance_distances(L);
W=diag(diag(L))-L;
W(W<1e-10)=0;
W=sparse(W);
[start_nodes,end_nodes,weights]=find(tril(W));
weights=max(0,weights);
Re=max(0,resistance_distances(sub2ind(size(resistance_distances),start_nodes,end_nodes)));
Pe=weights.*Re;
Pe=Pe/sum(Pe);

C0=1/30; % Rudelson, 1996 Random Vectors in the Isotropic Position (too hard to figure out actual C0)
C= 4*C0; % Rudelson and Vershynin, 2007, Thm. 3.1
q=round(9*C^2*N*log(N)/(epsilon^2));


results=gendist(Pe',q,1);
spin_counts=tabulate(results);
per_spin_weights=weights./(q*Pe);

counts=zeros(size(weights));
max_realized=max(results);
counts(1:max_realized)=spin_counts(:,2);
new_weights=counts.*per_spin_weights;
sparserW=sparse(start_nodes,end_nodes,new_weights,N,N);
sparserW=sparserW+sparserW';
sparserL=diag(sum(sparserW))-sparserW;
sparserD=diag(diag(sparserL));
sparserW=sparserD-sparserL;
sparserA=sign(sparserW);
end


