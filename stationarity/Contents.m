% GSPBOX - Stationarity   
%    
%  Vertex stationarity
%      gsp_experimental_psd - Experimental power density function
%      gsp_design_translates - Create a filterbank by uniformly translating a window
%      gsp_stationarity_cov - Covariance matrix from graph stationary data
%      gsp_stationarity_ratio - Assert the stationarity level of some data
%      gsp_estimate_psd - Estimation of the Power spectrum density
%      gsp_wiener_optimization - Solve wiener optimization problem
%      gsp_wiener_l2 - Solve wiener optimization problem with l2 fidelity term
%      gsp_wiener_inpainting - Solve wiener in-painting problem
%
%  For help, bug reports, suggestions etc. please send email to
%  gspbox 'dash' support 'at' groupes 'dot' epfl 'dot' ch
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/stationarity/Contents.php

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



