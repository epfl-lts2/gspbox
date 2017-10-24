% use a prior a*||W||_F^2 + b*||D||_F^2 instead of
%             a*||W||_F^2 + a*||D||_F^2 = a*||L||_F^2
%
% to be tried f.e. with demo_good_basis_index.m, demo_cumulative_energy.m 
%
% July 2016

function [W, stat] = gsp_learn_graph_l2_degrees_mine(Z, a, b, params)
%GSP_LEARN_GRAPH_L2_DEGREES Learn graph from pairwise distances using l2 prior on nodes degrees
%   Usage:  [W, stat] = gsp_learn_graph_l2_degrees(Z, a)
%           [W, stat] = gsp_learn_graph_l2_degrees(Z, a, params)
%
%   Inputs:
%         Z         : Matrix with (squared) pairwise distances of nodes
%         a         : ||L||_F^2 prior constant   (bigger a -> more dense W)
%         params    : Optional parameters
%
%   Outputs:
%         W         : Weighted adjacency matrix
%         stat      : Optional output statistics (adds small overhead)
%
%   'W = gsp_learn_graph_l2_degrees(Z, a, params)' computes a weighted
%   adjacency matrix $W$ from squared pairwise distances in $Z$, using the
%   smoothness assumption that $\text{trace}(X^TLX)$ is small, where $X$ is
%   the data (columns) changing smoothly from node to node on the graph and
%   $L = D-W$ is the combinatorial graph Laplacian.
%
%   Alternatively, Z can contain other types of distances and use the
%   smoothness assumption that
%
%   .. sum(sum(W .* Z))
%
%   .. math:: \sum_i\sum_j W_{ij}Z_{ij}
%
%   is small. 
%
%   The minimization problem solved is 
%
%   .. minimize_W sum(sum(W .* Z)) + a/2 * ||W||_F^2/2 
%   ..            + a/2 * ||sum(W)||_2^2 
%   ..            s.t. sum(sum(W)) == n
%
%   .. math:: \min_W \sum_i\sum_j W_{ij}Z_{ij} +\frac{\alpha}{2}\|W\|_F^2  +\frac{\alpha}{2}\|W1\|_2^2 s.t. \|W\|=n 
%
%   subject to $W$ being a valid weighted adjacency matrix (non-negative,
%   symmetric, with zero diagonal). Note that 
%
%   .. .||W||_F^2/2 + ||sum(W)||_2^2  = ||L||^2
%
%   .. math:: \|W\|_F^2 + \|W1\|_2^2 = \|L\|_F^2
%
%   The algorithm used is forward-backward-forward (FBF) based primal dual
%   optimization (see references).
%
%   Example:::
%
%         G = gsp_sensor(256);
%         f1 = @(x,y) sin((2-x-y).^2);
%         f2 = @(x,y) cos((x+y).^2);
%         f3 = @(x,y) (x-.5).^2 + (y-.5).^3 + x - y;
%         f4 = @(x,y) sin(3*((x-.5).^2+(y-.5).^2));
%         X = [f1(G.coords(:,1), G.coords(:,2)), f2(G.coords(:,1), G.coords(:,2)), f3(G.coords(:,1), G.coords(:,2)), f4(G.coords(:,1), G.coords(:,2))];
%         figure; subplot(2,2,1); gsp_plot_signal(G, X(:,1)); title('1st smooth signal');
%         subplot(2,2,2); gsp_plot_signal(G, X(:,2)); title('2nd smooth signal');
%         subplot(2,2,3); gsp_plot_signal(G, X(:,3)); title('3rd smooth signal');
%         subplot(2,2,4); gsp_plot_signal(G, X(:,4)); title('4th smooth signal');
%         Z = gsp_distanz(X').^2;
%         [W] = gsp_learn_graph_l2_degrees(Z*25, 1);
%         W(W<1e-5) = 0;
%         G2 = gsp_update_weights(G, W);
%         figure; gsp_plot_graph(G2); title('Graph with edges learned from above 4 signals');
%       
%
%   Additional parameters
%   ---------------------
%  
%   * *params.W_init*   : Initialization point. default: zeros(size(Z))
%   * *verbosity*       : Above 1 adds a small overhead. Default: 1
%   * *maxit*           : Maximum number of iterations. Default: 1000
%   * *tol*             : Tolerance for stopping criterion. Default: 1e-5
%   * *step_size*       : Step size from the interval (0,1). Default: 0.5
%
%   The stopping criterion is whether both relative primal and dual
%   distance between two iterations are below a given tolerance. 
%   
%   To set the step size use the following rule of thumb: Set it so that
%   relative change of primal and dual converge with similar rates (use
%   verbosity > 1).
%
%
%   See also: gsp_learn_graph_log_degrees gsp_distanz gsp_update_weights
% 
%   References: kalofolias2016learn komodakis2015playing kalofolias2017large
%

% Author: Vassilis Kalofolias
% Testing: gsp_test_learn_graph
% Date: June 2015


%% Default parameters
if nargin < 3
    params = struct;
end

if not(isfield(params, 'verbosity')),   params.verbosity = 1;   end
if not(isfield(params, 'maxit')),       params.maxit = 1000;      end
if not(isfield(params, 'tol')),         params.tol = 1e-5;      end
if not(isfield(params, 'step_size')),   params.step_size = .5;      end     % from (0, 1)


%% Fix parameter size and initialize
z = squareform(Z);
z = full(z(:));
l = length(z);                   % number of edges
% n(n-1)/2 = l => n = (1 + sqrt(1+8*l))/ 2
n = round((1 + sqrt(1+8*l))/ 2);    % number of nodes

if isfield(params, 'W_init')
    w = squareform(params.W_init)';
else
    % w_0 = exp(-6*dx2/mean(dx2));
    w = zeros(size(z));
end

%% Needed operators

% needed for gradient of h(w)
% L*w = sum(W)
[S, St] = sum_squareform(n);

% ones_l = ones(l,1);

sum_op = @(w) S*w;
sum_t_op = @(z) St*z;
% edges -> scalar
K_op = @(w) 2*sum(w);   % sum(sum(W)) = 2 * sum(w)
% scalar -> edges
% Kt_op = @(z) 2 * z * ones_l;
norm_K = 2 * sqrt(l);
%norm_L = l;         % do i need the square norm or the square root?

%% TODO: Rescaling??

%% Learn the graph
% min_{W>=0}     tr(X'*L*X) - gc * sum(log(sum(W))) + gp * norm(W-W0,'fro')^2, where L = diag(sum(W))-W
% min_W       I{W>=0} + W(:)'*Dx(:)  - gc * sum(log(sum(W))) + gp * norm(W-W0,'fro')^2
% min_W                f(W)          +       g(L_op(W))      +   h(W)

% put proximal of trace plus positivity together
f.eval = @(w) 2*w'*z;    % half should be counted
%f.eval = @(W) 0;
f.prox = @(w, c) max(0, w - 2*c*z);  % all change the same

% projection of sum of W on n
% g(w) = IND(sum(w) == n)
g.eval = @(z) 10 * abs(z-n);  % this could be 0. Measures violation of the condition sum(w) = n!!
g.prox = @(z, c) n;
% proximal of conjugate of g: z - c*g.prox(z/c, 1/c)
g_star_prox = @(z, c) z - c*n;

% the following is the squared frobenius norm of the graph Laplacian 
%      ... for b==a
h.eval = @(w) a * (2*norm(w)^2) + b*(norm(sum_op(w))^2);   % CAREFUL: w two times!
h.grad = @(w) 2 * (a*2*w + b*sum_t_op(sum_op(w)));
%  =  normest(4*a*speye(l) + 2*b*sum_t_op(sum_op(speye(l))))
h.beta = 4*a + 2*b*(2*(n-1));    % norm(ones(n) + a*eye(n)) = n+a

% if there is no quadratic term, it reduces to a linear program:
if a == 0
    %   X =     linprog(f,   A,   b,   Aeq,           beq,LB,UB,X0)
    [w, stat.fval, stat.exitflag, stat.output, stat.lambda] = linprog(z, [], [], ones(size(z))', 100*n, zeros(size(z)));
else
    %% My custom FBF based primal dual (see [Komodakis, Pesquet])
    % parameters mu, epsilon for convergence (see [1])
    mu = h.beta + norm_K;     %TODO: is it squared or not??
    epsilon = lin_map(0.0, [0, 1/(1+mu)], [0,1]);   % in (0, 1/(1+mu) )
    
    % INITIALIZATION
    % primal variable ALREADY INITIALIZED
    %w = params.w_init;
    % dual variable
    v_n = K_op(w);
    if nargout > 1 || params.verbosity > 1
        stat.f_eval = nan(params.maxit, 1);
        stat.g_eval = nan(params.maxit, 1);
        stat.h_eval = nan(params.maxit, 1);
        stat.fgh_eval = nan(params.maxit, 1);
        stat.pos_violation = nan(params.maxit, 1);
    end
    if params.verbosity > 1
        fprintf('Relative change of primal, dual variables, and objective fun\n');
    end
    tic
    gn = lin_map(params.step_size, [epsilon, (1-epsilon)/mu], [0,1]);              % in [epsilon, (1-epsilon)/mu]
    for i = 1:params.maxit
        %Y_n = w - gn * (h.grad(w) + Kt_op(v_n));
        Y_n = w - gn * (h.grad(w) + 2*v_n);
        y_n = v_n + gn * (K_op(w));
        P_n = f.prox(Y_n, gn);
        p_n = g_star_prox(y_n, gn); % = y_n - gn*g_prox(y_n/gn, 1/gn)
        %Q_n = P_n - gn * (h.grad(P_n) + Kt_op(p_n));
        Q_n = P_n - gn * (h.grad(P_n) + 2*p_n);
        q_n = p_n + gn * (K_op(P_n));
        
        if nargout > 1 || params.verbosity > 2
            stat.f_eval(i) = f.eval(w);
            stat.g_eval(i) = g.eval(K_op(w));
            stat.h_eval(i) = h.eval(w);
            stat.fgh_eval(i) = stat.f_eval(i) + stat.g_eval(i) + stat.h_eval(i);
            stat.pos_violation(i) = -sum(min(0,w));
        end
        rel_norm_primal = norm(- Y_n + Q_n, 'fro')/norm(w, 'fro');
        rel_norm_dual = norm(- y_n + q_n)/norm(v_n);
        
        if params.verbosity > 2
            fprintf('iter %4d: %6.4e   %6.4e   %6.4e    %6.4e   %6.3e \n', i, rel_norm_primal, rel_norm_dual, stat.fgh_eval(i), full(sum(w)*2), stat.h_eval(i));
        elseif params.verbosity > 1
            fprintf('iter %4d: %6.4e   %6.4e\n', i, rel_norm_primal, rel_norm_dual);
        end
        
        if i <= 1 && isfield(params, 'W_init') 
            w = squareform(params.W_init)';
        else
            w = w - Y_n + Q_n;
        end
        v_n = v_n - y_n + q_n;
        
        if rel_norm_primal < params.tol && rel_norm_dual < params.tol
            break
        end
    end
    stat.time = toc;
    if params.verbosity > 0
        fprintf('# iters: %4d. Rel primal: %6.4e Rel dual: %6.4e   %6.3e\n', i, rel_norm_primal, rel_norm_dual, f.eval(w) + g.eval(K_op(w)) + h.eval(w));
        fprintf('Time needed is %f seconds\n', stat.time);
    end
    
    % Use the following for testing:
    %     g.L = K_op;
    %     g.Lt = Kt_op;
    %     g.norm_L = norm_K;
    %     [w, info] = fbf_primal_dual(w, f, g, h, params);
    %[w, info] = fb_based_primal_dual(w, f, g, h, params);
    %%
    
    if params.verbosity > 3
        figure; plot(real([stat.f_eval, stat.g_eval, stat.h_eval])); hold all; plot(real(stat.fgh_eval), '.'); legend('f = w''*dx2', 'g = 1000 * abs(2sum(w)-n)', 'h = a*norm(L)^2', 'f+g+h');
        figure; plot(stat.pos_violation); title('sum of negative (invalid) values per iteration')
        figure; semilogy(max(0,-diff(real(stat.fgh_eval'))),'b.-'); hold on; semilogy(max(0,diff(real(stat.fgh_eval'))),'ro-'); title('|f(i)-f(i-1)|'); legend('going down','going up');
    end
end

W = squareform(w);

