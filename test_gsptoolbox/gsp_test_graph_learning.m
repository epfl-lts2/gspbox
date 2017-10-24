function n_errors = gsp_test_graph_learning


tol = 1e-5;

n_errors = test_graph_l2(tol) + test_graph_log(tol); %+ test_graph_LX(tol*100) ;
% n_errors = n_errors + test_graph_LX(tol*100);

end

function n_errors = test_graph_LX(tol)

n_errors = 0;

n = 20;
b = 1;

X = rand(20, 5);

params.verbosity = 0;
params.maxit = 10000;

W = gsp_learn_graph_LX_fro(X, b, params);

% zero diagonal
if not(all(diag(W) == 0))
    n_errors = n_errors + 1;
    warning('W has non zero diagonal');
end

% non-negative weights
if not(all(W(:) >= -tol))
    n_errors = n_errors + 1;
    warning('W has non-negative values');
end

% sum equal to n
if (norm(sum(W) - b) / n) > tol
    n_errors = n_errors + 1;
    warning('W has degrees smaller than the limit requested');
end

% W symmetric
if (norm(W - W', 'fro') / n) > tol
    n_errors = n_errors + 1;
    warning('W not symmetric');
end


end


function n_errors = test_graph_l2(tol)

n_errors = 0;

n = 20;

Z = gsp_distanz(randn(n)).^2;


params.maxit = 100000;
params.tol = 0.1 * tol;
W = gsp_learn_graph_l2_degrees(1*Z, 1, params);

% zero diagonal
if not(all(diag(W) == 0))
    n_errors = n_errors + 1;
end

% non-negative weights
if not(all(W(:) >= -tol))
    n_errors = n_errors + 1;
end

% sum equal to n
if abs(sum(W(:) - n) / n > tol)
    n_errors = n_errors + 1;
end
  
% % step-size correct
% params.step_size = .9;
% W2 = gsp_learn_graph_l2_degrees(10*Z, 10, params);
% if not(norm(W - W2, 'fro')/norm(W, 'fro') < 10*tol)
%     n_errors = n_errors + 1;
% end





end





function n_errors = test_graph_log(tol)

n_errors = 0;

n = 20;

Z = gsp_distanz(rand(n)).^2;

% We'll get a very sparse graph and check the degrees quality later
W = gsp_learn_graph_log_degrees(Z, 1, 1);

params.tol = tol * 1e-3;
params.maxit = 50000;
params.step_size = .9;
params.verbosity = 1;
W = gsp_learn_graph_log_degrees(Z, .031, .003, params);

% zero diagonal
if not(all(diag(W) == 0))
    n_errors = n_errors + 1;
end

% non-negative weights
if not(all(W(:) >= -tol))
    n_errors = n_errors + 1;
end

% positive degrees
if not(all(sum(W) >= tol * max(sum(W))))
    n_errors = n_errors + 1;
end


end





