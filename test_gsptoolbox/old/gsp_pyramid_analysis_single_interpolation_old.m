function [coarse_approximation,prediction_error]=gsp_pyramid_analysis_single_interpolation(signal,G,keep_inds,h_filter,param)
%GSP_PYRAMID_ANALYSIS_SINGLE_INTERPOLATION Compute a single level of the graph pyramid transform coefficients 
%   Usage:  [coarse_approximation,prediction_error]=gsp_pyramid_analysis_single(signal,G,keep_inds,h_filter,param);
%
%   Input parameters:
%         signal                    : Graph signal to analyze.
%         G                         : Graph structure on which the signal resides.
%         keep_inds                 : The indices of the vertices to keep when downsampling the graph and signal.
%         h_filter                  : The H filter in the Laplacian pyramid.
%         param                     : Contains optional additional parameters to be passed on to gsp_filter and gsp_interpolate (param.regularize_epsilon and param.use_exact for gsp_interpolate).
%   Output parameters:
%         coarse_approximation      : The next coarse approximation.
%         prediction_error          : The next prediction error.
%
%   'gsp__pyramid_analysis_single_interpolation(signal,G,keep_inds,h_filter,param)' computes 
%   a single level of the graph pyramid transform coefficients.
%
%   See also:  
%
%   Demos:  
% 
%   References: 

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:

% Compute the next coarse approximation and prediction error
h_filtered_signal=gsp_filter(G,h_filter,signal,param); 
coarse_approximation=h_filtered_signal(keep_inds);
prediction=gsp_interpolate_old(h_filtered_signal,G,keep_inds,param);
prediction_error=signal-prediction;