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
%   *param* is a structure containing optional arguments including
%   * *param.regularize_epsilon* : Interpolation parameter.
%
%   Please read the documentation of |gsp_filter_analysis| for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_analysis 
%             gsp_pyramid_synthesis gsp_pyramid_cell2coeff 
%
%   Demo: gsp_demo_pyramid
% 
%   References: shuman2013framework 

% Author: David I Shuman, Nathanael Perraudin
% Date: 26 November 2015
% Testing : test_pyramid

% Compute the next coarse approximation and prediction error
h_filtered_signal = gsp_filter(G, h_filter, signal, param); 
coarse_approximation = h_filtered_signal(keep_inds, :);
prediction = gsp_interpolate(G, coarse_approximation, keep_inds, param);
prediction_error = signal - prediction;