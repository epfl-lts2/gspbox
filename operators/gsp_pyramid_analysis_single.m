function [coarse_approximation,prediction_error]=gsp_pyramid_analysis_single(G,signal,keep_inds,h_filter,param)
%GSP_PYRAMID_ANALYSIS_SINGLE Compute a single level of the graph pyramid transform coefficients 
%   Usage:  [coarse_approximation,prediction_error]=gsp_pyramid_analysis_single(G,signal,keep_inds,h_filter,param);
%
%   Input parameters:
%         signal                    : Graph signal to analyze.
%         G                         : Graph structure on which the signal resides.
%         keep_inds                 : The indices of the vertices to keep when downsampling the graph and signal.
%         h_filter                  : The H filter in the pyramid transform.
%   Output parameters:
%         coarse_approximation      : The next coarse approximation.
%         prediction_error          : The next prediction error.
%
%   'gsp_pyramid_analysis_single(G,signal,keep_inds,h_filter,param)' computes 
%   a single level of the graph pyramid transform coefficients.
%   
%   param is a structure containing optional arguments including
%    param.regularize_epsilon : Interpolation parameter.
%
%   Please read the documentation of GSP_FILTER_ANALYSIS for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_analysis 
%             gsp_pyramid_synthesis gsp_pyramid_cell2coeff 
%
%   Demo: gsp_demo_pyramid
% 
%   References:
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_analysis_single.php

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

% Author: David I Shuman, Nathanael Perraudin
% Date: 26 November 2015
% Testing : test_pyramid

% Compute the next coarse approximation and prediction error
h_filtered_signal = gsp_filter(G, h_filter, signal, param); 
coarse_approximation = h_filtered_signal(keep_inds, :);
prediction = gsp_interpolate(G, coarse_approximation, keep_inds, param);
prediction_error = signal - prediction;
