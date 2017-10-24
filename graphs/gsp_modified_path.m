function [G]=gsp_modified_path(W)
%GSP_MODIFIED_PATH  Initialize a modified path graph
%   Usage:  G = gsp_modified_path(W);
%
%   Input parameters:
%         W     : Vector of weights.
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_modified_path(N)' initializes a path graph structure. The node are
%   connected like a path but with the weight given in W. The number of
%   node is the length of W + 1.
%
%   Example:
%
%          W = ones(15,1);
%          W(8) = 0.1;
%          G = gsp_modified_path(W);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_ring, gsp_path
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_modified_path.html

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

% Author: Nathanael Perraudin
% Date: 15 March 2014

N=length(W)+1;
G.N = N;
% Create weighted adjancency matrix
i_inds = [1:N-1,2:N];
j_inds = [2:N,1:N-1];
G.W = sparse(i_inds,j_inds,[W(:);W(:)],N,N);

% Create coordinates
d = zeros(N,1);
for ii=2:N
   d(ii) = d(ii-1)+log(1+1./W(ii-1))/log(2); 
end

G.coords = [d(:), zeros(N,1)];
G.limits = [-max(d)*0.05, max(d)*1.05, -1, 1];

G.type = 'modified path';

G = gsp_graph_default_parameters(G);

end

