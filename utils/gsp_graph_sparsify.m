function [ Gnew] = gsp_graph_sparsify(G,epsilon)
%GSP_GRAPH_SPARSIFY sparsify a graph using Spielman-Srivastava algorithm
%   Usage: Gnew = gsp_graph_sparsify(G,epsilon);
%
%   Input parameters
%       G       : Graph structure or Laplacian matrix
%       epsilon : Sparsification parameter
%
%   Ouput parameters:
%       Gnew    : New sparsified graph or new laplacian
%
%   This function sparsifies a graph using Spielman-Srivastava algorithm.
%   Note that epsilon should be between 1/sqrt(N) and 1.
%
%   Example:
%
%         epsilon = 0.5;
%         param.distribute = 1;
%         nnparam.k = 20;
%         param.nnparam = nnparam;
%         G = gsp_sensor(256,param);
%         G2 = gsp_graph_sparsify(G,epsilon);
%         figure(1)
%         gsp_plot_graph(G);
%         title('Original graph')
%         figure(2)
%         gsp_plot_graph(G2);
%         title('Sparsified graph')
%
%   References:
%     D. A. Spielman and N. Srivastava. Graph sparsification by effective
%     resistances. SIAM Journal on Computing, 40(6):1913--1926, 2011.
%     
%     M. Rudelson. Random vectors in the isotropic position. Journal of
%     Functional Analysis, 164(1):60--72, 1999.
%     
%     M. Rudelson and R. Vershynin. Sampling from large matrices: An approach
%     through geometric functional analysis. Journal of the ACM (JACM),
%     54(4):21, 2007.
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_graph_sparsify.php

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

% Author: David Shuman, Nathanael Perraudin
% Date  : 22 June 2014
% Testing: test_sparsify


% Test the input parameters
if isstruct(G)
    if isfield(G,'lap_type')
        if ~strcmp(G.lap_type,'combinatorial')
            error('Not implemented yet');
        end
    end
    L = G.L;
else
    L = G;
end 

N=size(L,1);

if N<3
    error('GSP_GRAPH_SPARSIFY: Cannot sparsify a graph with less than 3 nodes');
end

% Epsilon should be between 1/sqrt(N) and 1, with a larger epsilon leading to a sparser graph
if ( (epsilon <= 1/sqrt(N)) || (epsilon >1) )
    error('GSP_GRAPH_SPARSIFY: Epsilon out of required range');
end


% Compute resistance distances
resistance_distances = gsp_resistance_distances(L);

% Get the weight matrix, and check if the original graph is connected
if isstruct(G)
    W = G.W;
else
    W = diag(diag(L)) - L;
end
W(W<1e-10) = 0;
W = sparse(W);
original_connected=gsp_check_connectivity_undirected(W);
if (original_connected==0)
    warning('Original graph not connected before sparsification');
end

% Initialize the probability distribution that will be used to select edges in the sparsified graph
[start_nodes,end_nodes,weights]=find(tril(W));
weights=max(0,weights);
Re=max(0,resistance_distances(sub2ind(size(resistance_distances),start_nodes,end_nodes)));
Pe=weights.*Re;
Pe=Pe/sum(Pe);
    
max_tries=10; % maximum number of tries to get a connected graph

for i=1:max_tries
    
    % Set Q
    C0=1/30; % Rudelson, 1996 Random Vectors in the Isotropic Position (too hard to figure out actual C0)
    C= 4*C0; % Rudelson and Vershynin, 2007, Thm. 3.1
    q=round(9*C^2*N*log(N)/(epsilon^2));

    % Choose random edges in the new graph according the probability
    % distribution initialized above
    results=gendist(Pe',q,1);
    spin_counts=hist(results,1:length(Pe'));
    per_spin_weights=weights./(q*Pe);

    % Tally the new weights and form the new graph
    new_weights=spin_counts'.*per_spin_weights;
    
    sparserW=sparse(start_nodes,end_nodes,new_weights,N,N);
    sparserW=sparserW+sparserW';
    sparserL=diag(sum(sparserW))-sparserW;
    sparserA=sign(sparserW);

    % Check if new graph is connected. If not, reduce epsilon and try again
    new_graph_connected=gsp_check_connectivity_undirected(sparserA);
    if new_graph_connected
        break;
    elseif i==max_tries
        warning('Despite attempts to reduce epsilon, sparsified graph is disconnected');
    else
        epsilon=epsilon-(epsilon-1/sqrt(N))/2;
    end
end

if isstruct(G)
    if ~G.directed
        sparserW = (sparserW + sparserW')/2;
        sparserL = (sparserL + sparserL')/2;
    end
    Gnew=G;
    Gnew.W = sparserW;
    Gnew.L = sparserL;
    Gnew.A = sparserA;
    Gnew.d = diag(sparserL);
    Gnew.Ne = nnz(Gnew.W)/2;
else
    Gnew = sparse(sparserL);
end


end


