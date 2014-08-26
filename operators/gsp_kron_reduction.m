function Gnew = gsp_kron_reduction( G,ind )
%GSP_KRON_REDUCTION Compute the kron reduction
%   Usage: Gnew = gsp_kron_reduction( G,ind );
%          Wnew = gsp_kron_reduction( W,ind );
%
%   Input parameters:
%       G       : Graph structure or weight matrix
%       ind     : indices of the nodes to keep
%   Output parameters:
%       Gnew    : New graph structure or weight matrix       
%   
%   This function perform the Kron reduction of the weight matrix in the
%   graph G, with boundary nodes labeled by ind. This function will
%   create a new graph with a weight matrix Wnew that contain only boundary
%   nodes and is computed as the Schur complement of the original matrix
%   with respect to the selected indices.
%
%   If a matrix is given, then a martrix is returned
%
%   Example:
%
%           N = 64;
%           param.distribute = 1;
%           param.Nc = 5;
%           param.regular = 1;
%           G = gsp_sensor(N,param);
%           ind = 1:2:N;
%           Gnew = gsp_kron_reduction( G,ind );
%           figure;
%           subplot(121)
%           gsp_plot_graph(G);
%           title('Original graph');
%           subplot(122)
%           gsp_plot_graph(Gnew);
%           title('Kron reduction');
%       
%
%   References:
%     F. Dorfler and F. Bullo. Kron reduction of graphs with applications to
%     electrical networks. Circuits and Systems I: Regular Papers, IEEE
%     Transactions on, 60(1):150-163, 2013.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_kron_reduction.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

% Author: Nathanael Perraudin, Pierre Vandergheynst
% Date  : 23 July 2014
% Testing: test_operators


if isstruct(G)
    if ~strcmp(G.lap_type,'combinatorial')
        error('Not implemented');
    end
    if G.directed
        error('GSP_KRON_REDUCTION: this only work for undirected graphs.');
    end
    L = G.L;
else
    L = G;
end


N = size(L,1);
ind_comp = setdiff(1:N,ind);

L_red = L(ind, ind);
L_in_out = L(ind, ind_comp);
L_out_in = L(ind_comp, ind);
L_comp = L(ind_comp, ind_comp);

Lnew = L_red - L_in_out * (L_comp \ L_out_in);

% Make the laplacian symetric if it is almost symetric!
if sum(sum(abs(Lnew-transpose(Lnew))))<eps*sum(sum(abs(Lnew)))
    Lnew = (Lnew+transpose(Lnew))/2;
end

if isstruct(G)
    % Suppress the diagonal ? This is a good question?
    Wnew = diag(diag(Lnew)) - Lnew;
    Snew = diag(Lnew) - transpose(sum(Wnew));
    if norm(Snew)<eps(1000)
        Snew = 0;
    end
    Wnew = Wnew + diag(Snew);
    
    Gnew = gsp_copy_graph_attributes(G,0);
    Gnew.plotting.limits = G.plotting.limits;
    Gnew.coords = G.coords(ind,:);
    Gnew.W = Wnew;
    Gnew.type='Kron reduction';
    Gnew = gsp_graph_default_parameters(Gnew);

else
    Gnew = Lnew;
end

end


