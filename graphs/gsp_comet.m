function [G] = gsp_comet(N,k)
%GSP_COMET Initialize a comet graph
%   Usage:  G = gsp_comet(N,k);
%
%   Input parameters:
%         N     : Number of vertices. (default 32)
%         k     : Degree of center vertex. (default 12)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_comet(N,k)' initializes the commet graph. The commet graph is a
%   simple path graph with a star od degree k at its end.
%
%   Example:
%
%          G = gsp_comet(16,8);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%   See also: gsp_graph, gsp_ring, gsp_path
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_comet.html

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

% Author : David I Shuman, Nathanael Perraudin
% Date : 15 March 2015

if nargin < 2
    k = 12;
end

if nargin < 1
    N = 32;
end


% Create weighted adjancency matrix
i_inds=[ones(1,k),2:k+1,k+1:N-1,k+2:N];
j_inds=[2:k+1,ones(1,k),k+2:N,k+1:N-1];
G.W=sparse(i_inds,j_inds,ones(1,length(i_inds)),N,N);

% Create coordinates
G.coords=zeros(N,2);
inds=(2:k+1)';
G.coords(2:k+1,:)=[cos((inds-1)*2*pi/k),sin((inds-1)*2*pi/k)];
G.coords(k+2:end,1)=(2:N-k)';
G.plotting.limits=[-2,max(G.coords(:,1)),min(G.coords(:,2)),max(G.coords(:,2))];

G.type = 'comet';

G = gsp_graph_default_parameters(G);

end

