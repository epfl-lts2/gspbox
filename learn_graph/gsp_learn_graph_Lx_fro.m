function [W, stat] = gsp_learn_graph_Lx_fro(X, params)


% Author: Vassilis Kalofolias
% Date: March 2016
%
%   Url: http://lts2research.epfl.ch/gsp/doc/learn_graph/gsp_learn_graph_Lx_fro.php

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


%% Default parameters
if nargin < 2
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

if isfield(params, 'w_init')
    w = params.w_init(:);
else
    % w_0 = exp(-6*dx2/mean(dx2));
    w = zeros(size(z));
end

%% Needed operators

% needed for gradient of h(w)
% L*w = sum(W)
[S, St] = LG_create_S_matrix(n);

ones_l = ones(l,1);

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
f.eval = @(w) 0;    % half should be counted
%f.eval = @(W) 0;
f.prox = @(w, c) max(0, w);  % just project to positive

% projection of sum of W on n
% g(w) = IND(sum(w) == n)
g.eval = @(z) 1000 * abs(z-n);  % this could be 0. Measures violation of the condition sum(w) = n!!
g.prox = @(z, c) n;
% proximal of conjugate of g: z - c*g.prox(z/c, 1/c)
g_star_prox = @(z, c) z - c*n;

% the following is the squared frobenius norm of the graph Laplacian * X
% h.eval = @(w) a * (2*norm(w)^2 + norm(sum_op(w))^2);   % CAREFUL: w two times!
% h.grad = @(w) 2 * a * (2*w + sum_t_op(sum_op(w)));
% h.beta = 2 * a * (n+1);    % = norm(ones(n)+eye(n))
h.eval = @(w) a * (norm(w)^2 + norm(sum_op(w))^2);   % CAREFUL: w two times!
h.grad = @(w) 2 * a * (2*w + sum_t_op(sum_op(w)));
h.beta = 2 * a * (n+1);    % = norm(ones(n)+eye(n))

% if there is no quadratic term, it reduces to a linear program:
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
    Y_n = w - gn * (h.grad(w) + v_n);
    y_n = v_n + gn * (K_op(w));
    P_n = f.prox(Y_n, gn);
    p_n = g_star_prox(y_n, gn); % = y_n - gn*g_prox(y_n/gn, 1/gn)
    %Q_n = P_n - gn * (h.grad(P_n) + Kt_op(p_n));
    Q_n = P_n - gn * (h.grad(P_n) + p_n);
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
        fprintf('iter %3d: %6.4e   %6.4e   %6.3e    %6.3e   %6.3e \n', i, rel_norm_primal, rel_norm_dual, stat.fgh_eval(i), full(sum(w)*2), stat.h_eval(i));
    elseif params.verbosity > 1
        fprintf('iter %3d: %6.4e   %6.4e\n', i, rel_norm_primal, rel_norm_dual);
    end
    
    w = w - Y_n + Q_n;
    v_n = v_n - y_n + q_n;
    
    if rel_norm_primal < params.tol && rel_norm_dual < params.tol
        break
    end
end
stat.time = toc;
if params.verbosity > 0
    fprintf('# iters: %3d. Rel primal: %6.4e Rel dual: %6.4e   %6.3e\n', i, rel_norm_primal, rel_norm_dual, f.eval(w) + g.eval(K_op(w)) + h.eval(w));
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


W = squareform(w);


