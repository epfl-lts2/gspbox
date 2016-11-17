function [W, stat] = gsp_learn_graph_log_degrees(Z, a, b, params)
%GSP_LEARN_GRAPH_LOG_DEGREES Learn graph from pairwise distances using negative log prior on nodes degrees
%   Usage:  [W, stat] = gsp_learn_graph_log_degrees(Z, a, b)
%           [W, stat] = gsp_learn_graph_log_degrees(Z, a, b, params)
%
%   Inputs:
%         Z         : Matrix with (squared) pairwise distances of nodes
%         a         : Log prior constant  (bigger a -> bigger weights in W)
%         b         : W||_F^2 prior constant  (bigger b -> more dense W)
%         params    : Optional parameters
%
%   Outputs:
%         W         : Weighted adjacency matrix
%         stat      : Optional output statistics (adds small overhead)
%
%   
%   'W = gsp_learn_graph_log_degrees(Z, a, b, params)' computes a weighted
%   adjacency matrix W from squared pairwise distances in Z, using the
%   smoothness assumption that text{trace}(X^TLX) is small, where X is
%   the data (columns) changing smoothly from node to node on the graph and
%   L = D-W is the combinatorial graph Laplacian. See the paper of the
%   references for the theory behind the algorithm.
%
%   Alternatively, Z can contain other types of distances and use the
%   smoothness assumption that
%
%      sum(sum(W .* Z))
%
%   is small. 
%
%   The minimization problem solved is 
%
%      minimize_W sum(sum(W .* Z)) - a * sum(log(sum(W))) + b * ||W||_F^2/2
%
%   subject to W being a valid weighted adjacency matrix (non-negative,
%   symmetric, with zero diagonal).
%
%   The algorithm used is forward-backward-forward (FBF) based primal dual
%   optimization (see references).
%
%   Example:
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
%         % we can multiply the pairwise distances with a number to control sparsity
%         [W] = gsp_learn_graph_log_degrees(Z*25, 1, 1);
%         % clean up zeros
%         W(W<1e-5) = 0;
%         G2 = gsp_update_weights(G, W);
%         figure; gsp_plot_graph(G2); title('Graph with edges learned from above 4 signals');
%       
%
%   Additional parameters
%   ---------------------
%  
%    params.W_init   : Initialization point. default: zeros(size(Z))
%    verbosity       : Default = 1. Above 1 adds a small overhead
%    maxit           : Maximum number of iterations. Default: 1000
%    tol             : Tolerance for stopping criterion. Defaul: 1e-5
%    step_size       : Step size from the interval (0,1). Default: 0.5
%
%   The function is using the UNLocBoX functions sum_squareform and
%   squareform_sp. 
%   The stopping criterion is whether both relative primal and dual
%   distance between two iterations are below a given tolerance. 
%   
%   To set the step size use the following rule of thumb: Set it so that
%   relative change of primal and dual converge with similar rates (use
%   verbosity > 1).
%
%   See also: gsp_learn_graph_l2_degrees gsp_distanz gsp_update_weights  
% 
%   References:
%     V. Kalofolias. How to learn a graph from smooth signals. Technical
%     report, AISTATS 2016: proceedings at Journal of Machine Learning
%     Research (JMLR)., 2016.
%     
%     N. Komodakis and J.-C. Pesquet. Playing with duality: An overview of
%     recent primal? dual approaches for solving large-scale optimization
%     problems. Signal Processing Magazine, IEEE, 32(6):31--54, 2015.
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/learn_graph/gsp_learn_graph_log_degrees.php

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

% Author: Vassilis Kalofolias
% Testing: gsp_test_learn_graph
% Date: June 2015


%% Default parameters
if nargin < 4
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

w_0 = 0;

%% Needed operators
% S*w = sum(W)
[S, St] = sum_squareform(n);

% S: edges -> nodes
K_op = @(w) S*w;

% S': nodes -> edges
Kt_op = @(z) St*z;

% the next is an upper bound if params.fix_zeros
norm_K = sqrt(2*(n-1));


%% TODO: Rescaling??
% we want    h.beta == norm_K   (see definition of mu)
% we can multiply all terms by s = norm_K/2*b so new h.beta==2*b*s==norm_K


%% Learn the graph
% min_{W>=0}     tr(X'*L*X) - gc * sum(log(sum(W))) + gp * norm(W-W0,'fro')^2, where L = diag(sum(W))-W
% min_W       I{W>=0} + W(:)'*Dx(:)  - gc * sum(log(sum(W))) + gp * norm(W-W0,'fro')^2
% min_W                f(W)          +       g(L_op(W))      +   h(W)

% put proximal of trace plus positivity together
f.eval = @(w) 2*w'*z;    % half should be counted
%f.eval = @(W) 0;
f.prox = @(w, c) max(0, w - 2*c*z);  % all change the same

param_prox_log.verbose = params.verbosity - 3;
g.eval = @(z) -a * sum(log(z));
g.prox = @(z, c) prox_sum_log(z, c*a, param_prox_log);
% proximal of conjugate of g: z-c*g.prox(z/c, 1/c)
g_star_prox = @(z, c) z - c*a * prox_sum_log(z/(c*a), 1/(c*a), param_prox_log);

h.eval = @(w) b * norm(w - w_0,'fro')^2;
h.grad = @(w) 2 * b * (w - w_0);
h.beta = 2 * b;


%% My custom FBF based primal dual (see [1] = [Komodakis, Pesquet])
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
    Y_n = w - gn * (h.grad(w) + Kt_op(v_n));
    y_n = v_n + gn * (K_op(w));
    P_n = f.prox(Y_n, gn);
    p_n = g_star_prox(y_n, gn); % = y_n - gn*g_prox(y_n/gn, 1/gn)
    Q_n = P_n - gn * (h.grad(P_n) + Kt_op(p_n));
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

    if params.verbosity > 3
        fprintf('iter %4d: %6.4e   %6.4e   %6.3e', i, rel_norm_primal, rel_norm_dual, stat.fgh_eval(i));
    elseif params.verbosity > 2
        fprintf('iter %4d: %6.4e   %6.4e   %6.3e\n', i, rel_norm_primal, rel_norm_dual, stat.fgh_eval(i));
    elseif params.verbosity > 1
        fprintf('iter %4d: %6.4e   %6.4e\n', i, rel_norm_primal, rel_norm_dual);
    end
    
    % use a few iterations with same w to initialize the duals as well
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
    fprintf('# iters: %4d. Rel primal: %6.4e Rel dual: %6.4e  OBJ %6.3e\n', i, rel_norm_primal, rel_norm_dual, f.eval(w) + g.eval(K_op(w)) + h.eval(w));
    fprintf('Time needed is %f seconds\n', stat.time);
end

% use the following for testing:
% g.L = K_op;
% g.Lt = Kt_op;
% g.norm_L = norm_K;
% [w, info] = fbf_primal_dual(w, f, g, h, params);
% %[w, info] = fb_based_primal_dual(w, f, g, h, params);


%%

if params.verbosity > 3
    figure; plot(real([stat.f_eval, stat.g_eval, stat.h_eval])); hold all; plot(real(stat.fgh_eval), '.'); legend('f', 'g', 'h', 'f+g+h');
    figure; plot(stat.pos_violation); title('sum of negative (invalid) values per iteration')
    figure; semilogy(max(0,-diff(real(stat.fgh_eval'))),'b.-'); hold on; semilogy(max(0,diff(real(stat.fgh_eval'))),'ro-'); title('|f(i)-f(i-1)|'); legend('going down','going up');
end

W = squareform(w);



