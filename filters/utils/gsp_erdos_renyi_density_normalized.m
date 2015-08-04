function [pdf]=gsp_erdos_renyi_density_normalized(y,N,p)
%GSP_ERDOS_RENYI_DENSITY_NORMALIZED  The asymptotic density function of the normalized graph Laplacian
%eigenalues of an Erdos-Renyi random graph
%   Usage:  pdf=gsp_erdos_renyi_density_normalized(y,N,p);
%           
%   Input parameters:
%         y               : A vector of values at which the density will be evaluated.
%         N               : Number of vertices in the graph.
%         p               : Edge probability.
%   Output parameters:
%         pdf             : The values of the density function at the points y.
%
%   'gsp_erdos_renyi_density_normalized(y,N,p)' computes the asymptotic density
%   function (for large N) of the normalized graph Laplacian
%   eigenvalues of an Erdos-Renyi random graph with N vertices and edge
%   probability p. This is a semicircular density function.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/utils/gsp_erdos_renyi_density_normalized.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
pdf=sqrt((p*N)/(1-p))*(1/(2*pi))*real((y>1-2*sqrt((1-p)/(p*N))).*(y<1+2*sqrt((1-p)/(p*N))).*sqrt(4-p*N/(1-p)*(1-y).^2));

