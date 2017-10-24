function [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(Gs,signal,num_levels,param)
%GSP_PYRAMID_ANALYSIS Compute the graph pyramid transform coefficients 
%   Usage:  [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(Gs,signal,num_levels);
%           [coarse_approximations,prediction_errors]=gsp_pyramid_analysis(Gs,signal,num_levels,param);
%
%   Input parameters:
%         Gs                      : A multiresolution sequence of graph structures, including the idx parameters tracking the subsampling pattern.
%         signal                  : Graph signal to analyze.
%         num_levels              : Number of levels in the pyramid transform.
%   Output parameters:
%         coarse_approximations   : Cell array with the coarse approximations at each level.
%         prediction_errors       : Cell array with the prediction errors at each level.
%
%   'gsp_pyramid_analysis(Gs,signal,num_levels)' computes 
%   the graph pyramid transform coefficients of a signal $f$.
%   
%   *param* is a structure containing optional arguments including
%
%   * *param.regularize_epsilon* : Interpolation parameter.
%   * *param.h_filters* : A cell array of graph spectral filters. If just
%     one filter is included, it is used at every level of the pyramid. 
%     Default 
%
%     ..     h(x) = 0.5 / ( 0.5 + x)
%   
%     .. math:: h(x) = \frac{0.5}{0.5+x}
%
%   Please read the documentation of |gsp_filter_analysis| for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_synthesis 
%             gsp_pyramid_cell2coeff gsp_pyramid_analysis_single
%
%   Demo: gsp_demo_pyramid
% 
%   References: shuman2013framework 

% Author: David I Shuman, Nathanael Perraudin
% Date: 26 November 2015
% Testing : test_pyramid
  
% Read input parameters and check that inputs have the correct sizes
if nargin < 4
    param = struct;
end


if length(signal) ~= Gs{1}.N
    error('The signal to analyze should have the same dimension as the first graph');
end

if num_levels >= length(Gs)
    error('Not enough graphs provided to compute that many levels of the graph Laplacian pyramid');
end

if ~isfield(param,'h_filters')
    h_filters=cell(num_levels,1);
    for ii=1:num_levels
        h_filters{ii}=@(x) .5./(.5+x);
    end
elseif length(param.h_filters)==1
    h_filters=cell(num_levels,1);
    for ii=1:num_levels
        h_filters{ii}=param.h_filters;
    end
elseif length(param.h_filters)==num_levels
    h_filters=param.h_filters;
else
    error('param.h_filters should be a cell array of length 1 or num_levels');
end
        
% Compute the pyramid transform
coarse_approximations=cell(num_levels+1,1);
coarse_approximations{1}=signal;
prediction_errors = cell(num_levels,1);

for ii=1:num_levels
    [coarse_approximations{ii+1},prediction_errors{ii}] = ...
        gsp_pyramid_analysis_single( Gs{ii}, ...
        coarse_approximations{ii}, Gs{ii+1}.mr.idx, h_filters{ii}, param);
end

end