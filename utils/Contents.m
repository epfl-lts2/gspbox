% GSPBOX - Utils
%
%  Tests graphs
%    gsp_check_connectivity       - Check the connectivity
%    gsp_check_connectivity_undirected - Check the connectivity for undirected graph
%    gsp_isdirected               - Check if the graph is directed
%    gsp_check_weights            - Check the weight matrix
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
%   Graph operations
%    gsp_adj2vec                  - Precompute the gradient matrix
%    gsp_compute_fourier_basis    - Compute the graph Fourier basis
%    gsp_create_laplacian         - Create the laplacian from the weight matrix
%    gsp_estimate_lmax            - Estimate the maximum Laplacian eigenvalue
%    gsp_graph_sparsify           - Sparsify the graph
%    gsp_symmetrize               - Symmetrize a graph
%
%   Others
%    gsp_repmatline               - Variation of the function repmat
%    gsp_classic2graph_eig_order  - Mapping for the ordering of the ring and the DFT
%    gsp_reset_seed               - Set the seed for reproducibility
%    gsp_plotfig                  - Helping function to save figures
%    gsp_point2dcdf               - Points to discrete continuous density function
%    gsp_ddf2dcdf                 - Discrete to dicrete cumulative density function
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/Contents.php

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

% To be done
% - gsp_spectrum_cdf_approx
%  - gsp_tree_depths

