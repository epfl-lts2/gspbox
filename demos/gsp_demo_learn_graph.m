%GSP_DEMO_LEARN_GRAPH Demonstration of learning a graph from data
%
%   In this demo, we show how the graph learning can be used to learn a
%   graph from smoothly changing signals. The theory behind the algorithm
%   can be found in
%
%   [1] V. Kalofolias, How to learn a graph from smooth signals, AISTATS
%   2016.
%
%   [2] V. Kalofolias, N. Perraudin, Large Scale Graph Learning From Smooth Signals, ICLR
%   2019.
%
%   Suppose that we have some 2 dimensional smooth functions::
%
%             f1 = @(x,y) 20 * (-sin((2-x-y).^2)/2 + cos(y*3));
%             f2 = @(x,y) 30 * cos((x+y).^2);
%             f3 = @(x,y) 30 * ((x-.5).^2 + (y-.5).^3 + x - y);
%             f4 = @(x,y) 50 * sin(3*((x-.5).^2+(y-.5).^2));
%
%   and we have uniform samples as features, displayed below::
%
%             figure;
%             subplot(2,2,1); scatter(xc, yc, 700, X(:,1), '.'); 
%             title('1st smooth signal'); axis off; colorbar;
%             subplot(2,2,2); scatter(xc, yc, 700, X(:,2), '.');
%             title('2nd smooth signal'); axis off; colorbar;
%             subplot(2,2,3); scatter(xc, yc, 700, X(:,3), '.'); 
%             title('3rd smooth signal'); axis off; colorbar;
%             subplot(2,2,4); scatter(xc, yc, 700, X(:,4), '.'); 
%             title('4th smooth signal'); axis off; colorbar;
%
%   .. figure::
%
%      Different signals
%
%      
%
%   We can compute the pairwise distances of the features and learn a graph
%   using them::
%
%             Z1 = gsp_distanz(X(:, 1)').^2;
%             W1 = gsp_learn_graph_log_degrees(Z1, 1.5, 1, params);
%
%   The second parameter penalizes the formation of un-connected nodes, and
%   the third penalizes the formation of too strong weights. We then clean
%   any tiny edges (due to numerical error), to obtain a sparse weighted
%   adjacency matrix. We feed this to |gsp_graph| to create a graph with
%   the given coordinates and weights:: 
%             
%             W1(W1<1e-5) = 0;
%             G1 = gsp_graph(W1, [xc, yc]);
%
%   We can also update the weights of an already existing graph using
%   |gsp_update_weights|. If we learn the graphs of all four above
%   functions, we get quite different results::
%
%             figure;
%             subplot(2,2,1); gsp_plot_edges(G1, params_plot);
%             title('graph learned from 1st smooth signal');
%             subplot(2,2,2); gsp_plot_edges(G2, params_plot);
%             title('graph learned from 2nd smooth signal');
%             subplot(2,2,3); gsp_plot_edges(G3, params_plot);
%             title('graph learned from 3rd smooth signal');
%             subplot(2,2,4); gsp_plot_edges(G4, params_plot);
%             title('graph learned from 4th smooth signal');
%
%   .. figure::
%
%      Different graphs learned
%
%      
%
%   Note that the edges follow the level curves of the above functions.
%
%   If we use all four above smooth functions as features to learn the
%   graph::
%
%             Z = gsp_distanz(X').^2;
%             W = gsp_learn_graph_log_degrees(Z/500, 2, 1, params);
%
%   we get a result that has more local edges::
%
%             params_plot.show_edges = 1;
%             G.plotting.vertex_size = 5;
%             figure; gsp_plot_graph(G, params_plot);
%             title('Graph with edges learned from above 4 signals');
%             
%   .. figure::
%
%      Graph with edges learned from above 4 signals
%
%      
%
%   This is close to the graph that we would learn using the acutal
%   coordinates as features. So why does it work so well? We can see that
%   the pattern of the pairwise distances using these features is similar
%   to the one of the pairwise geometric distances between nodes::
%
%             figure;
%             subplot(1, 2, 1); 
%             imagesc(gsp_distanz(X')); 
%             title('Geometric pairwise distances between nodes');
%             subplot(1, 2, 2);
%             imagesc(gsp_distanz([xc, yc]'));
%             title('Pairwise distances computed from features');
%
%   .. figure::
%
%      Geometric pairwise distances between nodes
%
%      
%
%   .. figure::
%
%      Pairwise distances computed from features
%
%      
%
%   The functions available for learning a graph are
%   |gsp_learn_graph_log_degrees| and |gsp_learn_graph_l2_degrees|.
%
%
%   See also: gsp_learn_graph_l2_degrees gsp_demo_learn_graph_large
%
%   References: kalofolias2016learn kalofolias2017large


%% Initialization
clear
gsp_reset_seed(0);
method = 1;     %1 = logarithmic prior, 2 = L-2 prior
n = 100;

%% Sample space and define smooth signals
% we sample [0,1]^2 uniformly
coords = sortrows(rand(n, 2));
xc = coords(:, 1);
yc = coords(:, 2);

% these functions change smoothly on the 2D plane
f1 = @(x,y) 20 * (-sin((2-x-y).^2)/2 + cos(y*3));
f2 = @(x,y) 30 * cos((x+y).^2);
f3 = @(x,y) 30 * ((x-.5).^2 + (y-.5).^3 + x - y);
f4 = @(x,y) 50 * sin(3*((x-.5).^2+(y-.5).^2));
% f1 = @(x,y) exp((-sin((2-x-y).^2)/2 + cos(y*3)))*50;
% f2 = @(x,y) exp(cos((x+y).^2))*80;
% f3 = @(x,y) exp((x-.5).^2 + (y-.5).^3 + x - y)*80;
% f4 = @(x,y) exp(sin(3*((x-.5).^2+(y-.5).^2)))*100;

%% The features are the values of the above smooth functions on our samples
X = [f1(xc, yc), f2(xc, yc), f3(xc, yc), f4(xc, yc)];


%% Plot signals
figure; 
subplot(2,2,1); scatter(xc, yc, 700, X(:,1), '.'); title('1st smooth signal'); axis off; colorbar;
subplot(2,2,2); scatter(xc, yc, 700, X(:,2), '.'); title('2nd smooth signal'); axis off; colorbar;
subplot(2,2,3); scatter(xc, yc, 700, X(:,3), '.'); title('3rd smooth signal'); axis off; colorbar;
subplot(2,2,4); scatter(xc, yc, 700, X(:,4), '.'); title('4th smooth signal'); axis off; colorbar;

%%
params.maxit = 50000;
params.step_size = 0.1;
params.verbosity = 1;
params.tol = 1e-5;

if method == 1
    s = sqrt(2*(n-1))/2 / 3;
else
    s =  1/2/sqrt(2);
end
%%

Z1 = gsp_distanz(X(:, 1)').^2;
if method == 1
     W1 = gsp_learn_graph_log_degrees(s*Z1, s*1.5, s*1, params);
%    W1 = sqrt(1.5)*gsp_learn_graph_log_degrees(Z1/sqrt(1.5), 1, 1, params);
else
    W1 = gsp_learn_graph_l2_degrees(s*Z1, s*1.2, params);
end
W1(W1<1e-5) = 0;
G1 = gsp_graph(W1, [xc, yc]);

%%
Z2 = gsp_distanz(X(:, 2)').^2;
if method == 1
    W2 = gsp_learn_graph_log_degrees(s*Z2, s*3, s*1, params);
else
    W2 = gsp_learn_graph_l2_degrees(s*Z2, s*1, params);
end
W2(W2<1e-5) = 0;
G2 = gsp_graph(W2, [xc, yc]);

%%
Z3 = gsp_distanz(X(:, 3)').^2;
if method == 1
    W3 = gsp_learn_graph_log_degrees(s*Z3, s*5.5, s*1, params);
else
    W3 = gsp_learn_graph_l2_degrees(s*Z3, s*2, params);
end
W3(W3<1e-5) = 0;
G3 = gsp_graph(W3, [xc, yc]);

%%
Z4 = gsp_distanz(X(:, 4)').^2;
if method == 1
    W4 = gsp_learn_graph_log_degrees(s*Z4, s*1.9, s*1, params);
else
    W4 = gsp_learn_graph_l2_degrees(s*Z4, s*1, params);
end
W4(W4<1e-5) = 0;
G4 = gsp_graph(W4, [xc, yc]);

%%
Z = gsp_distanz(X').^2;
if method == 1
    W = gsp_learn_graph_log_degrees(s*Z/300, s*2, s*1, params);
else
    W = gsp_learn_graph_l2_degrees(s*Z/300, s*1, params);
end
W(W<1e-5) = 0;
G = gsp_graph(W / sum(W(:)) * n, [xc, yc]);

%%
% nnz(G0.W)
fprintf('Graph of signal 1: %d  edges\n', nnz(W1)/2);
fprintf('Graph of signal 2: %d  edges\n', nnz(W2)/2);
fprintf('Graph of signal 3: %d  edges\n', nnz(W3)/2);
fprintf('Graph of signal 4: %d  edges\n', nnz(W4)/2);
fprintf('Graph of all signals: %d  edges\n', nnz(W)/2);


%% Plot individual graphs from each signal
 params_plot.edge_size = 0;
% params_plot.show_edges = 0;
figure;
subplot(2,2,1); gsp_plot_edges(G1, params_plot); title('graph learned from 1st smooth signal'); ind = sum(G1.W>0)==0; hold on; scatter(xc(ind), yc(ind), 40, 'r');
subplot(2,2,2); gsp_plot_edges(G2, params_plot); title('graph learned from 2nd smooth signal'); ind = sum(G2.W>0)==0; hold on; scatter(xc(ind), yc(ind), 40, 'r');
subplot(2,2,3); gsp_plot_edges(G3, params_plot); title('graph learned from 3rd smooth signal'); ind = sum(G3.W>0)==0; hold on; scatter(xc(ind), yc(ind), 40, 'r');
subplot(2,2,4); gsp_plot_edges(G4, params_plot); title('graph learned from 4th smooth signal'); ind = sum(G4.W>0)==0; hold on; scatter(xc(ind), yc(ind), 40, 'r');

%% Plot graph from all signals
params_plot.edge_size = 1;
params_plot.show_edges = 1;
G.plotting.vertex_size = 5;
figure; gsp_plot_graph(G, params_plot); title('Graph with edges learned from above 4 signals');
hold on; ind = sum(G.W>0)==0; hold on; scatter(xc(ind), yc(ind), 40, 'r');


%% Why does it work so well?
figure; imagesc(gsp_distanz(X')); title('Geometric pairwise distances between nodes');
figure; imagesc(gsp_distanz([xc, yc]')); title('Pairwise distances computed from features');

% [dc2_sorted, ind_sorted] = sort(squareform(gsp_distanz([xc, yc]')), 'ascend');
% dx2 = squareform(gsp_distanz(X'));

% figure; plot(dc2_sorted/sum(dc2_sorted)); hold on; plot(dx2(ind_sorted)/sum(dx2));
% title('Sorted pairwise distances');
% legend('geometric distances', 'distances computed from features');

%% Other choices of smooth functions:
% f1 = @(x,y) (-sin((2-x-y).^2)/2 + cos(y*3))*50;
% f2 = @(x,y) cos((x+y).^2)*80;
% f3 = @(x,y) ((x-.5).^2 + (y-.5).^3 + x - y)*80;
% f4 = @(x,y) sin(3*((x-.5).^2+(y-.5).^2))*150;

% f1 = @(x,y) exp(2*((-sin((2-x-y).^2)/2 + cos(y*3))))*25;
% f2 = @(x,y) exp(2*(cos((x+y).^2)))*24;
% f3 = @(x,y) exp(2*((x-.5).^2 + (y-.5).^3 + x - y))*24;
% f4 = @(x,y) exp(2*(sin(3*((x-.5).^2+(y-.5).^2))))*30;

