%GSP_DEMO_GRAPH_EMBEDDING Demonstration file for graph embeddings
%
%   In this demo we perform a low-dimensional embedding on a specific graph
%   using three different algorithms and we plot the resulting embeddings
%   on a single plot. At first we create a 2-D sensor graph embedded in a
%   3-D space. We then compute 3 different embeddings for this graph:
%   Laplacian eigenmaps, LLE, and Isomap.
%
%   Laplacian eigenmaps
%   -------------------
%
%   This function uses the Laplacian of the graph and the diagonal degree
%   matrix to compute the eigenvalues and eigenvectors of generalized
%   eigenvector problem, L U =  D U Lambda. The coordinates of the
%   embedding are the eigenvectors that correspond to the bottom dim*
%   eigenvalues (ignoring always the zero eigenvalue).
%
%   Locally Linear Embedding (LLE) 
%   ------------------------------
%
%   This method first converts the weighed adjacency matrix to a distance
%   matrix. Then it computes another weight matrix A by minimizing each
%   reconstruction error epsilon = X -SUM_J A_J X_J ^2 where the
%   x_j are the neighboors of x. To do so we first compute the
%   following local covariance matrix:
%
%      C_{jk} = \frac{1}{2} (D_{:, j}+D_{i, :}-D_{ij}-D_{:,:})
%
%   where D_{ij} is the squared distance matrix,
%   
%      D_{:, j} = \frac{1}{n} \sum_i D_{ij},
%      D_{i,:} = \frac{1}{n} \sum_j D_{ij},
%      D_{:,:} = \frac{1}{n^2} \sum_i \sum_j D_{ij},
%
%   and where D is the n by n squared distance matrix (D = d^2).
%   One can observe that the local covariance matrix is simply a
%   Multi-Dimensional Scaling (MDS) of the squared distance matrix D.
%
%   We calculate the weight matrix A by solving the system of linear
%   equations sum_k C_{jk}a_k=1 and rescale the weights so that any
%   column of A sums up to 1. Finally the coordinates of the embedding are
%   the eigenvectors of the matrix M = (I - A)^T(I - A). More precisely
%   the coordinates of the embedding are the eigenvectors that correspond
%   to the dim+1 bottom eigenvalues of M (we always ignore the first
%   eigenvector since the first eigenvalue is always zero) leaving us with
%   a dim dimensional embedding.
%
%   Isomap
%   ------
%
%   This algorithm computes the embedding using the distance matrix d.
%   Firstly we compute the dijkstra distance of all possible nodes of the
%   graph. We store these values squared in a matrix D where D_{ij} is
%   the squared dijkstra distance from node i to node j. Continuing we
%   performe a MDS on the squared dijkstra matrix D. This way we compute
%   Matrix B as
%
%      B_{jk} = 0.5 * (D_{:, j}+D_{i,:}-D_{ij}- D_{:,:}). 
%
%   Finally the coordinates of the embedding are the eigenvectors that
%   corespond to the top dim eigenvalues. Finaly one can scale the
%   coordinates as X=Gamma Lambda ^{1/2} where Lambda is the
%   diagonal matrix of dim top eigenvalues of B and Gamma is the
%   matrix of the corresponding eigenvectors.
%
%   The signal on the graph is related to the coordinate information of the
%   original graph and therefore allows us to evaluate the resulting
%   embedding by looking at the smoothness of this signal on the graph. 
%
%   Figure 1: Resulting embeddings
%
%      
%
%   See also: gsp_nn_graph, gsp_compute_coordinates,
%   gsp_update_coordinates, gsp_plot_signal 
%
%   References:
%     L. K. Saul and S. T. Roweis. An introduction to locally linear
%     embedding. unpublished. Available at: http://www. cs. toronto. edu/~
%     roweis/lle/publications. html, 2000.
%     
%     M. Belkin and P. Niyogi. Laplacian eigenmaps and spectral techniques
%     for embedding and clustering. In NIPS, volume 14, pages 585--591, 2001.
%     
%     J. B. Tenenbaum, V. De Silva, and J. C. Langford. A global geometric
%     framework for nonlinear dimensionality reduction. Science,
%     290(5500):2319--2323, 2000.
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/demos/gsp_demo_graph_embedding.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

%% Initialization
clear
close all

% Call gsp_BOX
gsp_start
gsp_reset_seed(1)

% General parameters
n = 400;    % number of nodes
paramnn.k = 6; % number of nearest neighbors for the graph
dim = 2; % embedding dimentionality

%% Create a 2-d plane embeded in a 3-dimensional space

% Generate graph
alpha = rand(n, 2);
[Q, R] = qr(randn(3));
coords = alpha * Q(1:2,:);

G = gsp_nn_graph(coords,paramnn);
%% Compute new coordinates using different methods

new_c = gsp_compute_coordinates(G , dim, 'lle');
G2 = gsp_update_coordinates(G,new_c);

new_c2 = gsp_compute_coordinates(G , dim, 'laplacian_eigenmaps');
G3 = gsp_update_coordinates(G, new_c2);

new_c3 = gsp_compute_coordinates(G , dim, 'isomap');
G4 = gsp_update_coordinates(G, new_c3);

%% Plot the original coordinates on the new embeddings

figure
subplot(2,2,1);
gsp_plot_graph(G);
title('Original graph')

subplot(2,2,2);
title('LLE ')
gsp_plot_signal(G2,alpha(:,1));

subplot(2,2,3);
title('Laplacian eigenmaps')
gsp_plot_signal(G3,alpha(:,1));

subplot(2,2,4);
title('Isomap')
gsp_plot_signal(G4,alpha(:,1));

