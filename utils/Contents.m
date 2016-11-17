% GSPBOX - Utils
%
%  Tests graphs
%    gsp_check_connectivity       - Check the connectivity
%    gsp_check_connectivity_undirected - Check the connectivity for undirected graph
%    gsp_isdirected               - Check if the graph is directed
%    gsp_check_weights            - Check the weight matrix
%    gsp_check_fourier            - Check if the Fourier bais is computed
%    gsp_check_filtertype         - Check the filtertype
%
%   Norms
%    gsp_norm_tv                  - TV norm on graphs
%    gsp_norm_tik                 - Squared L2 norm of the gradient on graph
%    gsp_norm_l1_filterbank       - L1 norm of a signal after application of a filterbank
%    gsp_norm_l2_filterbank       - L2 norm of a signal after application of a filterbank
%    gsp_norm_tig                 - Norms of the atoms of a kernel
%
%   Distance
%    gsp_resistance_distances     - Compute the resistance distance
%    gsp_distanz                  - Compute all euclidean distances
%    gsp_nn_distanz               - Compute the nearest neighbors distances
%    gsp_rmse_mv                  - Compute all RMSE with missing values
%    gsp_hop_distanz              - Hop distances
%
%   Approximation
%    gsp_cheby_coeff              - Compute Chebysheff coefficients
%    gsp_cheby_op                 - Chebysheff polynomial approximation
%    gsp_cheby_eval               - Evaluate the Chebyshev polynomial
%    
%   Time-Vertex Approximation
%    gsp_jtv_cheby_coeff          - Compute Chebysheff coefficients of time-vertex filterbanks
%    gsp_jtv_cheby_op             - Chebysheff polynomial approximation of time-vertex filterbanks
%
%   Graph operations
%    gsp_adj2vec                  - Precompute the gradient matrix
%    gsp_vec2adj                  - Compute the weight marix from binary adjacency matrix and vector of weights
%    gsp_compute_fourier_basis    - Compute the graph Fourier basis
%    gsp_create_laplacian         - Create the laplacian from the weight matrix
%    gsp_incidence                - Create the incidence matrix
%    gsp_estimate_lmax            - Estimate the maximum Laplacian eigenvalue
%    gsp_fast_estimate_lk         - Estimate the k-th eigenvalue of the Laplacian
%    gsp_graph_sparsify           - Sparsify the graph
%    gsp_symmetrize               - Symmetrize a graph
%    gsp_assign_rand_direction    - Assign random direction to the edges of the graph
%
%   Others
%    gsp_jtv_fa                   - Create frequency axis
%    gsp_jtv_ta                   - Create time axis
%    gsp_decompose_asymmatrix     - Decompose a matrix in symmetric and asymmetric part
%    gsp_repmatline               - Variation of the function repmat
%    gsp_classic2graph_eig_order  - Mapping for the ordering of the ring and the DFT
%    gsp_reset_seed               - Set the seed for reproducibility
%    gsp_plotfig                  - Helping function to save figures
%    gsp_point2dcdf               - Points to discrete continuous density function
%    gsp_ddf2dcdf                 - Discrete to dicrete cumulative density function
%    gsp_good_graph_index         - Compute how well a graph fits a given data matrix
% 
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/Contents.php

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

% To be done
% - gsp_spectrum_cdf_approx
% - gsp_tree_depths

