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
%   ..   argmin   ||grad(x)||_1   s. t. ||Mx-b||_2 < epsilon
%
%   .. math:: arg \min_x  \|\nabla(x)\|_1 \text{ s. t. } \|Mx-b\|_2 \leq \epsilon
%
%   Where b represents the known measurements, M is an operator
%   representing the mask and $\epsilon$ is the radius of the l2 ball.
%
%   We set
%
%   * $f_1(x)=||\nabla x ||_1$
%     We define the prox of $f_1$ as:
%
%     .. prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||grad(z)||_1
%
%     .. math:: prox_{f1,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2 +  \gamma \| \nabla z \|_1
%
%   * $f_2$ is the indicator function of the set S define by $||Mx-b||_2 < \epsilon$
%     We define the prox of $f_2$ as
%
%     .. prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     .. math:: prox_{f2,\gamma} (z) = arg \min_{x} \frac{1}{2} \|x-z\|_2^2   + i_S(x) ,
%
%     with $i_S(x)$ is zero if x is in the set S and infinity otherwise.
%     This previous problem has an identical solution as:
%
%     .. argmin_{z} ||x - z||_2^2   s.t.  ||b - M z||_2 < epsilon
%
%     .. math:: arg \min_{z} \|x - z\|_2^2   \hspace{1cm} such \hspace{0.25cm} that \hspace{1cm} \|Mz-b\|_2 \leq \epsilon
%
%     It is simply a projection on the B2-ball.
%
%   Results
%   -------
%
%   .. figure::
%
%      Original signal on graph
%
%      This figure shows the original signal on graph.
%
%   .. figure::
%
%      Depleted signal on graph
%
%      This figure shows the signal on graph after the application of the
%      mask and addition of noise. Half of the vertices are set to 0.
%
%   .. figure::
%
%      Reconstructed signal on graph usign TV
%
%      This figure shows the reconstructed signal thanks to the algorithm.
%
%   Comparison with Tikhonov regularization
%   ---------------------------------------
%
%   We can also use the Tikhonov regularizer that will promote smoothness.
%   In this case, we solve:
%
%   ..   argmin   ||grad(x)||_2^2   s. t. ||Mx-b||_2 < epsilon
%
%   .. math:: arg \min_x \tau \|\nabla(x)\|_2^2 \text{ s. t. } \|Mx-b\|_2 \leq \epsilon
%
%   The result is presented in the following figure:
%
%   .. figure::
%
%      Reconstructed signal on graph using Tikhonov
%
%      This figure shows the reconstructed signal thanks to the algorithm.
%


% Author: Nathanael Perraudin
% Date: 4th March 2014


%% Initialisation
clear;
close all;
% Load toolbox
init_unlocbox();

%% Important parameters
% size of the graph for the demo
N = 250;
% probability of having a label on a vertex.
p = 0.1;
verbose = 1;    % verbosity level
sigma = 0.0;





%% Create a graph
%  paramgraph.distribute = 1;
% G = gsp_sensor(N, paramgraph);
G = gsp_sensor(N);
G = gsp_community(N);

% for a path with irregular weights and distances
% G = gsp_path(N);
% G = gsp_update_coordinates(G, G.coords + [rand(G.N, 1), zeros(G.N, 1)]);
% G = gsp_update_weights(G, G.W .* exp(-gsp_distanz(G.coords').^2));

G = gsp_adj2vec(G);
G = gsp_estimate_lmax(G);
G = gsp_compute_fourier_basis(G);

% x_0 = (1 + sign(G.U(:,4))) / 2;
x_0 = round((1 + sign(G.U(:,2))) * 2+(1 + sign(G.U(:,3))) +(1 + sign(G.U(:,4))) / 2);
% x_0 = round(linspace(0, 5, N))';
x_0 = G.info.node_com;

%% Set the labels
% create the mask
% ind_obs = boolean(full(sparse(1:4:G.N, 1, 1, G.N, 1)));
ind_obs = rand(G.N, 1) < p;
ind_unobs = not(ind_obs);

%applying the Mask to the data
x_obs = ind_obs.*(x_0+sigma*randn(G.N,1));


%% Method 1: Graph TV
% setting different parameter for the simulation
param_solver.verbose = verbose;  % display parameter
param_solver.tol = 1e-5;
param_solver.maxit = 4000;

tic;
sol_tv = gsp_regression_tv(G, ind_obs, x_obs, 0, param_solver);
toc

%% Method 2: Tikhonov
tic;
G2 = G;
G2.L = G.L^3;
%sol_tik = gsp_regression_tik(G, ind_obs, x_obs, 0, param_solver);
sol_tik = gsp_regression_tik(G2, ind_obs, x_obs, 0, param_solver);
toc

%% Method 3: My idea: inpaint using graph TV on one variable that is l-2
%% away from another one that is smooth
%%
% x, y are primal variables: z is the dual
%      min_x,y,z ||grad(x_all)||_1 + a/2||x-y||^2 + b/2 ||grad(y_all)||^2
%
% grad(x_all) = grad([x;x_obs]) = A*x + A_obs * x_obs = A*x - b
%
% only rows of gradient corresponding to unknown values used:
%      min_x,y,z ||Ax - b||_1  + a/2||x-y||^2 + b/2 ||Ay - b||^2
% Get rid of operator A by using dual variable:
%      min_x,y,z ||z - b||_1   + a/2||x-y||^2 + b/2 ||Ay - b||^2
% Think of all primal variables as one vector w = [x;y] containing ONLY
% the elements on the unlabeled nodes

% choose parameters for solving the problem:
alpha = .1;
beta = 1;

% compute the Diff matrix of the graph
G = gsp_adj2vec(G);

% operators based on the Diff
A_obs = G.Diff(:, ind_obs);
A_unobs = G.Diff(:, not(ind_obs));
n_unobs = nnz(not(ind_obs));
zeros_x = zeros(n_unobs, 1);
zeros_y = zeros_x;

% graph TV term is the l-1 on the dual:
b = - A_obs * x_obs(ind_obs);
f1_params.y = b;
f1.eval = @(z) norm(z - b, 1);
f1.prox = @(z, gamma) prox_l1(z, gamma, f1_params);
f1.L = @(w) A_unobs * w(1:n_unobs);
f1.Lt = @(z) [A_unobs' * z; zeros_y];
f1.norm_L = normest(A_unobs);
%f1.norm_L = sqrt(G.lmax);

% use prox for the second function
% w = [x;y], w0 = [x0;y0]
% prox_g(x) = ((1+g)x + gy)/(1+2g)
f2.prox = @(w, gamma) [ ((1+alpha*gamma) * w(1:n_unobs) + ...
    alpha*gamma     * w(n_unobs+1:end))/(1+2*alpha*gamma);...
    (alpha*gamma     * w(1:n_unobs) + ...
    (1+alpha*gamma) * w(n_unobs+1:end))/(1+2*alpha*gamma)];
f2.eval = @(w) alpha/2 * norm(w(1:n_unobs) - w(n_unobs+1:end))^2;

% use gradient for the last function
AtA = A_unobs' * A_unobs;
Atb = A_unobs' * b;
f3.grad = @(w) [zeros_x; beta * (AtA*w(n_unobs+1:end) - Atb)];
f3.eval = @(w) beta/2 * norm(A_unobs * w(n_unobs+1:end) - b)^2;
f3.beta = beta * normest(AtA);

%% solve the problem
param_solver.verbose = 1;
param_solver.algo = 'FBF_PRIMAL_DUAL';
param_solver.normalized_timestep = 0.99;
tic;
[sol_mine, info_sol] = solvep([sol_tv(not(ind_obs));  sol_tik(not(ind_obs))], {f1,f2,f3}, param_solver);
toc

sol_mine_x = x_obs;
sol_mine_x(ind_unobs) = sol_mine(1:n_unobs);
sol_mine_y = x_obs;
sol_mine_y(ind_unobs) = sol_mine(n_unobs+1:end);

%% Method 4: L-1 of Lx ||Lx||_1
tic;
param_solver = rmfield(param_solver, 'algo');
sol_Lx_l1 = gsp_regression_Lx_l1(G, ind_obs, x_obs, 0, param_solver);
toc


%% Compute the errors
ans_tv = round(sol_tv);
ans_tik = round(sol_tik);
ans_mine = round(sol_mine_x);
ans_Lxl1 = round(sol_Lx_l1);

err_tv      = nnz(ans_tv ~= x_0) / G.N
err_tik     = nnz(ans_tik ~= x_0) / G.N
err_mine    = nnz(ans_mine ~= x_0) / G.N
err_Lxl1    = nnz(ans_Lxl1 ~= x_0) / G.N

%% Print the result
paramplot.show_edges = 1;
val_lims = lin_map([-.05, 1.05], [min([x_0; sol_tv; sol_tik]), max([x_0; sol_tv; sol_tik])], [0, 1]);


%% Plot the original graph
if not(strcmp(G.type, 'path'))
    figure(1)
    paramplot.vertex_highlight = ind_obs;
    gsp_plot_signal(G, x_0, paramplot)
    caxis(val_lims);
    title('Original signal: highlighted known values')
    
    
    % % Let show depleted graph
    % figure(2)
    % gsp_plot_signal(G,depleted_graph_value,paramplot)
    % caxis([-1 1])
    % title('Measurement')
    % Let show the reconstructed graph
    figure(3)
    paramplot.vertex_highlight = (ans_tv ~= x_0);
    gsp_plot_signal(G, ans_tv, paramplot)
    caxis(val_lims)
    title('TV solution: highlighted mistakes')
    
    % Let show the reconstructed graph
    figure(4)
    paramplot.vertex_highlight = (ans_tik ~= x_0);
    gsp_plot_signal(G, ans_tik, paramplot)
    caxis(val_lims)
    title('Tikhonov solution: highlighted mistakes')
    
    % Let show the reconstructed graph
    figure(5)
    paramplot.vertex_highlight = (ans_mine ~= x_0);
    gsp_plot_signal(G, ans_mine, paramplot)
    caxis(val_lims)
    title('My solution: highlighted mistakes')
    
    % Let show the reconstructed graph
    figure(6)
    paramplot.vertex_highlight = (ans_Lxl1 ~= x_0);
    gsp_plot_signal(G, ans_Lxl1, paramplot)
    caxis(val_lims)
    title('||Lx||_1 solution: highlighted mistakes')
    
    %%
else strcmp(G.type, 'path')
    figure; 
    plot(G.coords(:, 1), sol_tv)
    hold on;
    plot(G.coords(:, 1), sol_tik)
    plot(G.coords(:, 1), x_0);
    plot(G.coords(:, 1), sol_Lx_l1)
    plot(G.coords(ind_obs, 1), x_0(ind_obs), 'o');
    legend('TV', 'Tikhonov', '||Lx||_1', 'true');
end


figure; plot(sol_mine_x)
hold on; plot(sol_mine_y)


