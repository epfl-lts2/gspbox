% GSPBOX - Operators
%
%  Localisation
%    gsp_translate      -  Translation operator
%    gsp_modulate       -  Modulation operator
%    gsp_localize       -  Localize a kernel
%
%  Differential
%    gsp_grad_mat       -  Compute the gradient sparse matrix
%    gsp_grad           -  Compute the gradient of a signal
%    gsp_div            -  Compute the divergence of a signal
%    gsp_adj2vec        -  Prepare the graph for the gradient computation
%
%  Transforms
%    gsp_gft            -  Graph Fourier transform
%    gsp_igft           -  Inverse graph Fourier transform
%    gsp_gwft           -  Windowed graph Fourier transform
%    gsp_ngwft          -  Normalized windowed graph Fourier transform
%
%  Pyramid - Reduction
%    gsp_kron_reduction - Kron reduction
%    gsp_kron_pyramid   - Compute the pyramide using kron reduction
%    gsp_pyramid_analysis - Analysis operator for graph pyramid
%    gsp_pyramid_synthesis - Sythesis operator for graph pyramid
%    gsp_pyramid_cell2coeff - Transform cell coefficient into a vector 
%    gsp_tree_multiresolution - Compute a tree reduction of a graph
%   
%
%   Utils
%    gsp_create_laplacian -  Create of change the laplacian of a graph
%    gsp_compute_fourier_basis - Compute the Fourier basis of a graph
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/index.php

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



