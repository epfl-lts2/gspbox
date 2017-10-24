function [reconstruction]=gsp_pyramid_synthesis_old(coarsest_approximation,prediction_errors,Gs,varargin)
%GSP_PYRAMID_SYNTHESIS Synthesizes a signal from its graph pyramid transform coefficients 
%   Usage:  reconstruction=gsp_pyramid_synthesis(coarsest_approximation,prediction_errors,Gs);
%           reconstruction=gsp_pyramid_synthesis(coarsest_approximation,prediction_errors,Gs,param);
%
%   Input parameters:
%         coarsest_approximation    : The coarsest approximation of the original signal.
%         prediction_errors         : Cell array with the prediction errors at each level.
%         Gs                        : A multiresolution sequence of graph structures, including the idx parameters tracking the subsampling pattern.
%   Output parameters:
%         reconstruction            : The synthesized signal.
%   Additional parameters:
%         param.use_exact           : To use exact graph spectral filtering instead of the Chebyshev approximation.
%         param.order               : Degree of the Chebyshev approximation (default=30).
%         param.least_squares       : Set to 1 to use the least squares synthesis (default=0) 
%         param.h_filters           : The filters used in the analysis operator. These are required for least squares synthesis, but not for the direct synthesis method
%
%   'gsp_pyramid_synthesis(coarsest_approximation,prediction_errors,Gs)' 
%   synthesizes a signal from its graph pyramid transform coefficients. 
%
%   See also:  
%
%   Demos:  
% 
%   References: 

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  

% Read input parameters 
if nargin>3
    param=varargin{1};
else
    param=0;
end

num_levels=length(prediction_errors);

% Compute the pyramid transform
current_coarse_approximation=coarsest_approximation;

for i=num_levels:-1:1
    if isfield(param, 'least_squares')
        if param.least_squares
            if ~isfield(param, 'h_filters')
                error('h-filter not provided');
            else
                param.h_filter=param.h_filters{i};
            end
        end
    end
    current_coarse_approximation=gsp_pyramid_synthesis_single_interpolation_old(current_coarse_approximation,prediction_errors{i},Gs{i},Gs{i+1}.idx,param);
end
reconstruction=current_coarse_approximation;


