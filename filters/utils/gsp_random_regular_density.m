function [pdf]=gsp_random_regular_density(y,r)
%GSP_RANDOM_REGULAR_DENSITY  The asymptotic density function of the graph Laplacian
%eigenalues of a random regular graph
%   Usage:  pdf=gsp_random_regular_density(y,r);
%           
%   Input parameters:
%         y               : A vector of values at which the density will be evaluated.
%         r               : The degree of every vertex in the graph.
%   Output parameters:
%         pdf             : The values of the density function at the points y.
%
%   'gsp_random_regular_density(y,r)' computes the asymptotic density
%   function (for large N) of the combinatorial graph Laplacian
%   eigenvalues of a random regular graph with each vertex having degree r.
%   This density function is given by McKay's Law.
%
%   See also:  
%
%   Demos:  
%
%   Requires: 
% 
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/utils/gsp_random_regular_density.php

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

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
pdf=real((y>r-2*sqrt(r-1)).*(y<r+2*sqrt(r-1)).*(r*sqrt(4*(r-1)-(r-y).^2)./(2*pi*(r^2-(r-y).^2))));
