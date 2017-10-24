%GSP_DEMO_GRAPH_TV Reconstruction of missing sample on a graph using TV
%
%   In this demo, we try to reconstruct missing sample of a piece-wise
%   smooth signal on a graph. To do so, we will minimize the well-known TV
%   norm defined on the graph.
%
%   For this example, you need the unlocbox. You can download it here:
%   http://unlocbox.sourceforge.net/download
%
%   We express the recovery problem as a convex optimization problem of the
%   following form:
%
%        argmin   ||grad(x)||_1   s. t. Mx = y
%
%   Where b represents the known measurements, M is an operator
%   representing the mask and epsilon is the radius of the l2 ball.
%
%   We set 
%
%    f_1(x)=||nabla x _1
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||grad(z)||_1
%
%    f_2 is the indicator function of the set S define by Mx = y
%     We define the prox of f_2 as 
%
%        prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     with i_S(x) is zero if x is in the set S and infinity otherwise.
%     This previous problem has an identical solution as:
%
%        argmin_{x} ||x - z||_2^2   s.t.  Mx = y
%
%     It is simply a projection on the B2-ball.
%
%   Results
%   -------
%
%   Figure 1: Original signal on graph
%
%      This figure shows the original signal on graph.
%
%   Figure 2: Depleted signal on graph
%
%      This figure shows the signal on graph after the application of the
%      mask and addition of noise. Half of the vertices are set to 0.
%
%   Figure 3: Reconstructed signal on graph usign TV
%
%      This figure shows the reconstructed signal thanks to the algorithm.
%
%   Comparison with Tikhonov regularization
%   ---------------------------------------
%
%   We can also use the Tikhonov regularizer that will promote smoothness.
%   In this case, we solve:
%   
%        argmin   ||grad(x)||_2^2   s. t. ||Mx-b||_2 < epsilon
%
%   The result is presented in the following figure:
%
%   Figure 4: Reconstructed signal on graph using Tikhonov
%
%      This figure shows the reconstructed signal thanks to the algorithm.
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/demos/gsp_demo_graph_tv.html

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
% Date: 4th March 2014


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 1;    % verbosity level
sigma = 0.0;

N = 256; % size of the graph for the demo





%% Create a random sensor graph

paramgraph.distribute = 1;
G = gsp_sensor(N,paramgraph);

G = gsp_adj2vec(G);
G = gsp_estimate_lmax(G);
G = gsp_compute_fourier_basis(G);

graph_value = sign(G.U(:,4));


%%
p = 0.6; %probability of having no label on a vertex.
%create the mask
M = rand(G.N,1);
M = M>p;


%applying the Mask to the data
depleted_graph_value = M.*(graph_value+sigma*randn(G.N,1));

sol = gsp_regression_tv(G,M,depleted_graph_value,0);
sol2 = gsp_regression_tik(G,M,depleted_graph_value,0);

%% Print the result
paramplot.show_edges = 1;

% Let show the original graph
figure(1)
gsp_plot_signal(G,graph_value,paramplot)
caxis([-1 1])
title('Original signal')


% Let show depleted graph
figure(2)
gsp_plot_signal(G,depleted_graph_value,paramplot)
caxis([-1 1])
title('Measurement')


% Let show the reconstructed graph
figure(3)
gsp_plot_signal(G,sol,paramplot)
caxis([-1 1])
title('Solution of the algorithm: TV')

% Let show the reconstructed graph
figure(4)
gsp_plot_signal(G,sol2,paramplot)
caxis([-1 1])
title('Solution of the algorithm: Tikhonov')

