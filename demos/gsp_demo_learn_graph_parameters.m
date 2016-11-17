n = 200;

%%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/demos/gsp_demo_learn_graph_parameters.php

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
X = 5 * rand(2, n);
Z = gsp_distanz(X).^2;

%% first way to control sparsity (Kalofolias 2016)
% this is easy to use when we solve an optimization problem that alternates
% between graph learning and solving w.r.t. X with fixed graph

alpha = .6;         % bigger weights
beta = .5;         % smaller big weights (more dense)
W1 = gsp_learn_graph_log_degrees(Z, alpha, beta);

%% second way to control sparsity (Kalofolias, Vandergheynst 2017)
% this is easy to use for controling sparsity. For a sparser graph just use
% a bigger theta. Also in general converges faster with the default
% stepsize. 
theta = 1/sqrt(alpha*beta);     % controls sparsity
delta = sqrt(alpha/beta);       % controls magnitude
W2 = delta * gsp_learn_graph_log_degrees(theta * Z, 1, 1);

%max(W2(:))
%max(W1(:))

%% the two ways are equivalent:
% Clean numerical error
W1(W1<1e-4)=0;
W2(W2<1e-4)=0;
fprintf('Relative error between W1 and W2: %e\n', norm(W1-W2, 'fro')/norm(W1, 'fro'));

%% Plot graph
% figure; spy(W1); title('Learned adjacency matrix');
% figure; spy(W2)
figure; gsp_plot_graph(gsp_graph(W1, X')); title('Learned graph');

