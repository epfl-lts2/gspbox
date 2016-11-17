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
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_classic2graph_eig_order.php

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


