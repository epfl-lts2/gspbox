n = 200;

%%
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
