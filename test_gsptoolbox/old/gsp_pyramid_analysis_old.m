function [coarse_approximations,prediction_errors]=gsp_pyramid_analysis_old(signal,Gs,num_levels,varargin)
%GSP_PYRAMID_ANALYSIS Compute the graph pyramid transform coefficients 
%   Usage:  [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(signal,Gs,num_levels);
%           [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(signal,Gs,num_levels,param);
%
%   Input parameters:
%         signal                    : Graph signal to analyze.
%         Gs                        : A multiresolution sequence of graph structures, including the idx parameters tracking the subsampling pattern.
%         num_levels                : Number of levels in the pyramid transform.
%   Output parameters:
%         coarse_approximations     : Cell array with the coarse approximations at each level.
%         prediction_errors         : Cell array with the prediction errors at each level.
%   Additional parameters:
%         param.use_exact           : To use exact graph spectral filtering instead of the Chebyshev approximation.
%         param.order               : Degree of the Chebyshev approximation (default=30).
%         param.regularize_epsilon  : Interpolation parameter.
%         param.h_filters           : A cell array of graph spectral filters. If just one filter is included, it is used at every level of the pyramid.
%
%   'gsp__pyramid_analysis(signal,Gs,num_levels)' computes 
%   the graph pyramid transform coefficients of a signal $f$.
%
%   See also:  
%
%   Demos:  
% 
%   References: 

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
% Read input parameters and check that inputs have the correct sizes
if nargin>3
    param=varargin{1};
else
    param=0;
end

if length(signal) ~= Gs{1}.N
    error('The signal to analyze should have the same dimension as the first graph');
end

if num_levels >= length(Gs)
    error('Not enough graphs provided to compute that many levels of the graph Laplacian pyramid');
end

if ~isfield(param,'h_filters')
    h_filters=cell(num_levels,1);
    for i=1:num_levels
        h_filters{i}=@(x) .5./(.5+x);
    end
elseif length(param.h_filters)==1
    h_filters=cell(num_levels,1);
    for i=1:num_levels
        h_filters{i}=param.h_filters;
    end
elseif length(param.h_filters)==num_levels
    h_filters=param.h_filters;
else
    error('param.h_filters should be a cell array of length 1 or num_levels');
end
        
% Compute the pyramid transform
coarse_approximations=cell(num_levels+1,1);
coarse_approximations{1}=signal;
prediction_errors=cell(num_levels,1);

for i=1:num_levels
    [coarse_approximations{i+1},prediction_errors{i}]=gsp_pyramid_analysis_single_interpolation_old(coarse_approximations{i},Gs{i},Gs{i+1}.idx,h_filters{i},param);
end
