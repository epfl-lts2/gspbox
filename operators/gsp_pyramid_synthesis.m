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
%   param is a structure containing optional arguments including
%
%    param.regularize_epsilon : Interpolation parameter.
%    param.least_squares : Set to 1 to use the least squares synthesis 
%     (default=0)  
%    param.use_landweber : Set to 1 to use the Landweber iteration 
%     approximation in the least squares synthesis.
%    param.landweber_its : Number of iterations in the Landweber 
%     approximation for least squares synthesis.
%    param.landweber_tau : Parameter for the Landweber iteration.
%    param.h_filters : A cell array of graph spectral filters. If just
%     one filter is included, it is used at every level of the pyramid. 
%     These filters are required for least squares synthesis, but not for 
%     the direct synthesis method 
%     Default 
%
%            h(x) = 0.5 / ( 0.5 + x)
%   
%   Please read the documentation of GSP_FILTER_ANALYSIS for other
%   optional arguments.
%
%   See also: gsp_graph_multiresolution gsp_pyramid_analysis 
%             gsp_pyramid_analysis_single gsp_pyramid_cell2coeff
%             gsp_pyramid_synthesis_single
%
%   Demo: gsp_demo_pyramid
% 
%   References:
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_synthesis.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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


