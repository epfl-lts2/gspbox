function [reconstruction,coarse_approximations]=gsp_pyramid_synthesis(Gs,coarsest_approximation,prediction_errors,param)
%GSP_PYRAMID_SYNTHESIS Synthesizes a signal from its graph pyramid transform coefficients 
%   Usage:  reconstruction=gsp_pyramid_synthesis(Gs,coarsest_approximation,prediction_errors);
%           reconstruction=gsp_pyramid_synthesis(Gs,coarsest_approximation,prediction_errors,param);
%
%   Input parameters:
%         Gs                      : A multiresolution sequence of graph structures, including the idx parameters tracking the subsampling pattern.
%         coarsest_approximation  : The coarsest approximation of the original signal.
%         prediction_errors       : Cell array with the prediction errors at each level.
%   Output parameters:
%         reconstruction          : The synthesized signal.
%         coarse_approximations   : Sequence of coarse approximations
%
%   'gsp_pyramid_synthesis(Gs,coarsest_approximation,prediction_errors)' 
%   synthesizes a signal from its graph pyramid transform coefficients. 
%  
%   *param* is a structure containing optional arguments including
%
%   * *param.regularize_epsilon* : Interpolation parameter.
%   * *param.least_squares* : Set to 1 to use the least squares synthesis 
%     (default=0)  
%   * *param.use_landweber* : Set to 1 to use the Landweber iteration 
%     approximation in the least squares synthesis.
%   * *param.landweber_its* : Number of iterations in the Landweber 
%     approximation for least squares synthesis.
%   * *param.landweber_tau* : Parameter for the Landweber iteration.
%   * *param.h_filters* : A cell array of graph spectral filters. If just
%     one filter is included, it is used at every level of the pyramid. 
%     These filters are required for least squares synthesis, but not for 
%     the direct synthesis method 
%     Default 
%
%     ..     h(x) = 0.5 / ( 0.5 + x)
%   
%     .. math:: h(x) = \frac{0.5}{0.5+x}
%
%   Please read the documentation of |gsp_filter_analysis| for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_analysis 
%             gsp_pyramid_analysis_single gsp_pyramid_cell2coeff
%             gsp_pyramid_synthesis_single
%
%   Demo: gsp_demo_pyramid
% 
%   References: shuman2013framework 

%   Author : David I Shuman, Nathanael Perraudin.
%   Date   : 26 November 2015
%   Testing: test_pyramid
  

% Read input parameters 
if nargin<4
    param=struct;
end
if ~isfield(param, 'least_squares') param.least_squares=0; end

num_levels=length(prediction_errors);

if param.least_squares
    if ~isfield(param, 'h_filters')
        error('h-filter not provided');
    end
end

if ~isfield(param,'h_filters')
    param.h_filters=cell(num_levels,1);
    for i=1:num_levels
        param.h_filters{i}=@(x) .5./(.5+x);
    end
elseif length(param.h_filters)==1
    if iscell(param.h_filters)
        filter=param.h_filters{1};
    else
        filter=param.h_filters;
    end
    param.h_filters=cell(num_levels,1);
    for i=1:num_levels
        param.h_filters{i}=filter;
    end
elseif length(param.h_filters)~=num_levels
    error('param.h_filters should be a cell array of length 1 or num_levels');
end


% Compute the inverse pyramid transform 
coarse_approximations = cell(num_levels+1,1);
coarse_approximations{num_levels+1} = coarsest_approximation;

for i=num_levels:-1:1
    param.h_filter=param.h_filters{i};
    coarse_approximations{i}=gsp_pyramid_synthesis_single(Gs{i},coarse_approximations{i+1},prediction_errors{i},Gs{i+1}.mr.idx,param);
end
reconstruction=coarse_approximations{1};

