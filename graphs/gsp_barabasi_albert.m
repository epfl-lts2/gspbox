function [ G ] = gsp_barabasi_albert( N, m0, m, param )
%GSP_ERDOS_RENYI Create a random Barabasi Albert graph
%   Usage:  G = gsp_barabasi_albert( N, m0, m, param );
%           G = gsp_barabasi_albert( N, m0 );
%           G = gsp_barabasi_albert( N );
%
%   Input parameters:
%         N     : Total number of nodes
%         m0    : Initial number of nodes
%         m     : Number of edges for each node addition. m <= m0.
%         param : Structure of optional parameter
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_barabasi_albert( N,p,param )' initializes create a graph structure
%   containing the weighted adjacency matrix (G.W), for a Barabasi Albert
%   graph. All edge weights are equal to 1. 
% 
%   The Barabasi Albert graph is constructed by connecting nodes randomly.
%   From a set of m0 initial nodes, one node is added a each step as well
%   as m edges from this new node to the already existing one following the
%   preferential attachment principle. This means that the probability of
%   creating a new edge from the new node to any other node is proportional
%   to its current degree.
%
%   param a Matlab structure containing the following fields:
% 
%      param.directed : defines if graph should be directed
%       By default, it is 0.
%
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%
%   Example:
%
%        G = gsp_barabasi_albert(100, 2, 2)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_barabasi_albert.php

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

% Author : Lionel Martin
 

% Optional input arguments
if nargin < 2, m0 = 1; end
if nargin < 3, m = 1; end
if nargin < 4, param = struct; end

if ~ isfield(param, 'directed'), param.directed = 0 ; end
if ~ isfield(param, 'verbose'), param.verbose = 1 ; end


% Check if the user didn't put silly parameter
if m0 < 1 || m < 1 || m0 < m
    error('The model requires 1 <= m <= m0.');
end

if N < m0
    error('N must be larger than m0.');
end

% Create the graph structure
G = struct;
G.graph_type = 'barabasi_albert';
G.type = 'barabasi_albert';

% Create the graph
G.W = sparse(N, N);

for id=m0+1:N
    distr = sum(G.W, 2);
    distr = distr + [ones(1, id-1), zeros(1, N-id+1)]';

    vidx = [];
    while 1
        is = randsample(N, m, true, distr);
        is = unique(is)

        distr(is) = 0;
        vidx = [vidx; is];
        size(vidx)

        if length(vidx) >= m
           for elem = 1:m
               G.W(vidx(elem), id) = 1;
           end
           break;
        end
    end
end

if ~param.directed
   G.W = G.W + G.W';
end

G = gsp_graph_default_parameters(G);

end

