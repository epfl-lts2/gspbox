function [G] = gsp_spiral(N,k,param)
%GSP_SPIRAL  Initialize a spriral graph
%   Usage:  G = gsp_spiral();
%           G = gsp_spiral(N);
%           G = gsp_spiral(N,k);
%           G = gsp_spiral(N, k, param);
%
%   Input parameters:
%         N     : Number of vertices. (default 200)
%         k     : Number of turns (3)
%         param : Structure of optional parameters
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_spiral(N)' initializes a graph structure for the sprial graph.
%   
%   The optional parameters are:
%    param.noise  : Noise level (Default 5)
%    param.start : Starting angle in degree (Default 90)
%    param.k : Approximate number of connections per nodes (Default 4)
%
%   Example:
%
%          G = gsp_spiral(64);
%          param.show_edges = 1;
%          gsp_plot_graph(G,param);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_spiral.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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

% Author : Nathanael Perraudin
% Date   : 11 December 2015

if nargin < 1
   N = 200; 
end

if nargin < 2
    k = 3;
end

if nargin<3
    param = struct;
end

if ~isfield(param,'noise'), param.noise = 5; end
if ~isfield(param,'start'), param.start =  90; end
if ~isfield(param,'K'), param.k =  4; end

deg2rad = (2*pi)/360;
start = param.start * deg2rad;

degrees = k*360;

n = start + sqrt(sort(rand(N,1))) * degrees * deg2rad;      
coords = [cos(n).*n+rand(N,1)*param.noise, -sin(n).*n+rand(N,1)*param.noise];

paramnn.k = 6*param.k;
Gt = gsp_nn_graph(coords,paramnn);
% keeps only the good connections
W = Gt.W;
W( logical( triu(ones(N),param.k) + tril(ones(N),-param.k ) ) )= 0;
 
G.W = W;

% Create coordinates
G.coords = coords;

G.type = 'spiral';

G = gsp_graph_default_parameters(G);

end

