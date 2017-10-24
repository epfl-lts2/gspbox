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
%   function (for large $N$) of the combinatorial graph Laplacian
%   eigenvalues of a random regular graph with each vertex having degree $r$.
%   This density function is given by McKay's Law.
%
%

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
pdf=real((y>r-2*sqrt(r-1)).*(y<r+2*sqrt(r-1)).*(r*sqrt(4*(r-1)-(r-y).^2)./(2*pi*(r^2-(r-y).^2))));
