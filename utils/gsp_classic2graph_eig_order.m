function [ v ] = gsp_classic2graph_eig_order( N )
%GSP_CLASSIC2GRAPH_EIG_ORDER Compute the graph eigen value ordering 
%   Usage: v = gsp_classic2graph_eig_order(N)
%       
%   Input parameters
%       N   : size of the graph
%   Output parameters
%       v   : vector of indexes
%
%   This function make the link between the DFT and the ring graph. It
%   returns the graph eigenvector ordering with respect of the DFT
%   ordering.
%

% Author: Nathanael Perraudin
% Date: 17 March 2014


if mod(N,2)
    % odd
    v = zeros(N,1);
    v(1) = 1;
    for ii = 1:((N-1)/2)
        v(2*ii) = ii+1;
        v(2*ii+1) = N-ii+1;
    end
else
    % even
    v = zeros(N,1);
    v(1) = 1;
    for ii = 1:(N/2-1)
        v(2*ii) = ii+1;
        v(2*ii+1) = N+1-ii;
    end
    v(N) = N/2+1;
end


end

