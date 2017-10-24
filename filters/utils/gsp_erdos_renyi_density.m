function [pdf]=gsp_erdos_renyi_density(y,N,p)
%GSP_ERDOS_RENYI_DENSITY  The asymptotic density function of the graph Laplacian
%eigenalues of an Erdos-Renyi random graph
%   Usage:  pdf=gsp_erdos_renyi_density(y,N,p);
%           
%   Input parameters:
%         y               : A vector of values at which the density will be evaluated.
%         N               : Number of vertices in the graph.
%         p               : Edge probability.
%   Output parameters:
%         pdf             : The values of the density function at the points y.
%
%   'gsp_erdos_renyi_density(y,N,p)' computes the asymptotic density
%   function (for large N) of the combinatorial graph Laplacian
%   eigenvalues of an Erdos-Renyi random graph with N vertices and edge
%   probability p. This is the density function associated with the free
%   additive convolution of a normal distribution with a semicircular
%   distribution.
%
%   See also:
%
%   Requires: gsp_free_conv_norm_semi
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/utils/gsp_erdos_renyi_density.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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
  
pdf=sqrt(1/((1-p)*N*p))*gsp_free_conv_norm_semi((y-p*N)/(sqrt(p*(1-p)*N))); 

