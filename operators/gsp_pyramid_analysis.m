function [ca,pe]=gsp_pyramid_analysis(Gs,f,param)
%GSP_PYRAMID_ANALYSIS Compute the graph pyramid transform coefficients 
%   Usage:  [ca,pe]=gsp_pyramid_analysis(Gs, f);
%           [ca,pe]=gsp_pyramid_analysis(Gs, f, param);
%
%   Input parameters:
%         Gs      : A multiresolution sequence of graph structures.
%         f       : Graph signal to analyze.
%         param   : Structure of optional parameters
%   Output parameters:
%         ca      : Cell array with the coarse approximation at each level
%         pe      : Cell array with the prediction errors at each level
%
%   'gsp_pyramid_analysis' computes the graph pyramid transform
%   coefficients of the signal f for the pyramid structure in Gs.
%
%   See also: gsp_kron_pyramid gsp_pyramid_synthesis gsp_pyramid_cell2coeff
%
%   Demo: gsp_demo_pyramid
% 
%   References:
%     I. Pesenson. Variational splines and paley-wiener spaces on
%     combinatorial graphs. Constructive Approximation, 29(1):1-21, 2009.
%     
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_analysis.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

% Author: Nathanael Perraudin
% Date: 5 August 2014
% Testing : test_pyramid



  

  
  
if nargin < 3
    param = struct;
end

if length(f) ~= Gs{1}.N
    error('The signal to analyze should have the same dimension as the first graph');
end


Nl = length(Gs)-1;



% Initisalization
ca = cell(Nl+1,1);
pe = cell(Nl+1,1);

ca{1} = f(:);
pe{Nl+1} = zeros(Gs{Nl+1}.N,1);

for ii = 2 : (Nl+1)
    % Low pass the signal
    s_low  = gsp_filter_analysis(Gs{ii-1},Gs{ii}.pyramid.filter,ca{ii-1},param);
    % Keep only the coefficient on the selected nodes
    ca{ii} = s_low(Gs{ii}.pyramid.ind);
    % Compute prediction
    s_pred = gsp_interpolate(Gs{ii-1},Gs{ii},ca{ii},param);
    % Compute errors
    pe{ii-1} = ca{ii-1} - s_pred;
end
    

end


