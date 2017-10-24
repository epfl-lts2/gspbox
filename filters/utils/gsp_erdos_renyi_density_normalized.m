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
%   function (for large $N$) of the normalized graph Laplacian
%   eigenvalues of an Erdos-Renyi random graph with $N$ vertices and edge
%   probability $p$. This is a semicircular density function.
%

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
pdf=sqrt((p*N)/(1-p))*(1/(2*pi))*real((y>1-2*sqrt((1-p)/(p*N))).*(y<1+2*sqrt((1-p)/(p*N))).*sqrt(4-p*N/(1-p)*(1-y).^2));
