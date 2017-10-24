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
%   function (for large $N$) of the combinatorial graph Laplacian
%   eigenvalues of an Erdos-Renyi random graph with $N$ vertices and edge
%   probability $p$. This is the density function associated with the free
%   additive convolution of a normal distribution with a semicircular
%   distribution.
%
%   See also:
%
%   Requires: gsp_free_conv_norm_semi
%

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
pdf=sqrt(1/((1-p)*N*p))*gsp_free_conv_norm_semi((y-p*N)/(sqrt(p*(1-p)*N))); 
