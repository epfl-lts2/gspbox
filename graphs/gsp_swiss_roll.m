function [G]=gsp_swiss_roll(N,s,thresh,rand_state)
%GSP_SWISS_ROLL Initialize a swiss roll graph
%   Usage:  G = gsp_swiss_roll(N,s,thresh,rand_state);
%
%   Input parameters:
%         N          : Number of vertices.
%         s          : sigma ( default: sqrt(2/N))
%         thresh     : threshold (default: 1e-6)
%         rand_state : rand seed (default: 45)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_create_swiss_roll(N,s,thresh,rand_state)' initializes a graph
%   structure containing the swiss roll graph
%
%   Example:
%
%          G = gsp_swiss_roll(200);
%          gsp_plot_graph(G);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_swiss_roll.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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
% Test: test_graphs

if nargin<1
    N = 200;
end

if nargin<2
   s = sqrt(2/N); 
end

if nargin<3
   thresh = 1e-6;
end

if nargin<4
   rand_state = 45; 
end


dataparams=struct('n',N,'dataset',-1','noise',0,'state',rand_state); % state was 0 before - 45 now
r=create_synthetic_dataset(dataparams);
G.coords=rescale_center(r.x)'; % r.x is a 3 x N matrix with the rows equal to the X, Y, and, Z components
G.plotting.limits = [-1,1,-1,1,-1,1];
dist=gsp_distanz(G.coords');
G.W=exp(-dist.^2/(2*s^2));
G.W = G.W-diag(diag(G.W));
G.W(G.W<thresh)=0;

G = gsp_graph_default_parameters(G);

end



