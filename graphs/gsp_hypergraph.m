function [G] = gsp_hypergraph(N,E, w, coords, limits)
%GSP_HYPERGRAPH  Initialize a hypergraph from a set of edges and weights
%   Usage:  G = gsp_hypergraph(W );
%           G = gsp_hypergraph(N,E);
%           G = gsp_hypergraph(N,E, w);
%           G = gsp_hypergraph(N,E, w, coords);
%           G = gsp_hypergraph(N, E, w, coords, limits);  
%
%   Input parameters:
%         N     : Number of nodes
%         E     : Set of edges (cell array)
%         w     : weights of the edges (default all ones)
%         coords: Coordonates of the points (optional)
%         limits: limits for the coordonates (optional)
%   Output parameters:
%         G     : Graph structure.
%
%   Example:
%
%         N = 100;
%         Nf = 2;
%         k = 4;
%         x = rand(N,Nf);
%         paramnn.k = k;
%         [indx, indy, d] = gsp_nn_distanz(x',x',paramnn);
%         sigma = mean(d)^2;
%         wt = exp(-d.^2/sigma);
%         E = cell(N,1);
%         w = zeros(N,1);
%         for ii = 1:N
%             edge = indx((1:k)+(ii-1)*k);
%             E{ii} = edge;
%             w(ii) = sum(wt(edge));
%         end
% 
%         G = gsp_hypergraph(N,E,w)
%
%   See also: gsp_nn_hypergraph gsp_graph
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_hypergraph.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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
% Date: 21  October 2015

if nargin<3
    w = ones(numel(E));
end

if nargin > 3
    G.coords = coords;
end
if nargin > 4
    G.plotting.limits = limits;
end
G.N = N;
G.Ne = numel(E);
G.W = sparse(N,G.Ne);
G.E = E;
for ii = 1:G.Ne
    % Here we use W for HW...
    G.W(G.E{ii},ii) = sqrt(w(ii));
end



G.type = 'hypergraph from edges';
G.directed = 0;
G.hypergraph = 1;

G = gsp_graph_default_parameters(G);

end

