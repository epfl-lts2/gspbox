% GSPBOX - Operators
%
%  Localisation
%    gsp_localize       -  Localize a kernel
%    gsp_modulate       -  Generalized modulation operator
%
%  Differential
%    gsp_grad_mat       -  Compute the gradient sparse matrix
%    gsp_grad           -  Compute the gradient of a signal
%    gsp_div            -  Compute the divergence of a signal
%
%  Transforms
%    gsp_gft            -  Graph Fourier transform
%    gsp_igft           -  Inverse graph Fourier transform
%    gsp_gwft           -  Windowed graph Fourier transform
%    gsp_ngwft          -  Normalized windowed graph Fourier transform
%
%  Time-Vertex Transforms
%    gsp_jft            -  Joint Time-Vertex Fourier transform
%    gsp_ijft           -  Inverse Joint Time-Vertex Fourier transform
%    gsp_tft            -  Time-Vertex Time-Fourier transform
%    gsp_itft           -  Time-Vertex Inverse Time-Fourier transform
%
%  Pyramid - Reduction
%    gsp_kron_reduce    - Kron reduction
%    gsp_graph_multiresolution - Compute a multiresolution of graphs
%    gsp_pyramid_analysis - Analysis operator for graph pyramid
%    gsp_pyramid_analysis_single - Compute a single level of the graph pyramid transform coefficients
%    gsp_pyramid_synthesis - Sythesis operator for graph pyramid
%    gsp_pyramid_synthesis_single -Synthesize a single level of the graph pyramid transform 
%    gsp_pyramid_cell2coeff - Keep only the necessary coefficients
%    gsp_interpolate    - Interpolate a signal
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
% see also: prox
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/Contents.php

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


