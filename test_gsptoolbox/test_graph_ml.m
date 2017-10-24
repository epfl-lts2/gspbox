function errors = test_graph_ml()

errors = 0;

errors = errors + test_gsp_regression_tik;
errors = errors + test_knn_classify_graph();
errors = errors + compare_flann_matlab();

%graph_ml_compare_all();

end


function errors = test_gsp_regression_tik()

errors = 0;

N = 50;
gsp_reset_seed(0)
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

x = randn(N,1);

y = gsp_filter(G,@(x) 1./(1+x),x);

M = rand(N,1) < 0.5;

param.direct = 1;
sol2 = gsp_regression_tik(G ,M, y , 0, param );

param.direct = 0;
param.maxit = 1000;
param.tol = 1e-13;
sol1 = gsp_regression_tik(G ,M, y , 0, param );



if norm(sol1-sol2)/norm(sol1) < 3e-5
    fprintf('GSP_PROX_TIK: exact ok\n');
else
    warning('GSP_PROX_TIK: exact Pas OK');
    badness = norm(sol1-sol2)/norm(sol1) 
    errors = errors+1;
end


M = ones(N,1);

tau = 5;
g = @(x) 1./(1+tau*x);
sol2 = gsp_filter(G,g,y);
param.maxit = 1000;
param.tol = 1e-13;
sol1 = gsp_regression_tik(G ,M, y , tau, param );


if norm(sol1-sol2)/norm(sol1) < 1e-5
    fprintf('GSP_PROX_TIK: filter ok\n');
else
    warning('GSP_PROX_TIK: filter Pas OK');
    badness = norm(sol1-sol2)/norm(sol1) 
    errors = errors+1;
end


% % G = gsp_sensor(N);
% % G = gsp_compute_fourier_basis(G);
% 
% x = randn(N,1);
% 
% y = gsp_filter(G,@(x) 1./(1+x),x);
% 
% M = rand(N,1) < 0.5;
% 
% param.direct = 1;
% param.exact = 1;
% sol2 = gsp_regression_tik(G ,M, y , 0, param );
% 
% param.exact = 0;
% param.order = 100;
% param.tol = 1e-13;
% sol1 = gsp_regression_tik(G ,M, y , 0, param );
% 
% 
% 
% if norm(sol1-sol2)/norm(sol1) < 1e-5
%     fprintf('GSP_PROX_TIK: exact 2 ok\n');
% else
%     warning('GSP_PROX_TIK: exact 2 Pas OK');
%     badness = norm(sol1-sol2)/norm(sol1) 
%     errors = errors+1;
% end
%%
    N = 100;
    G = gsp_sensor(N);
    G = gsp_estimate_lmax(G);
    M = rand(N,1)>0.5;
    tau = 2;
    
%    Mop =@(x) bsxfun(@times,M,x);
    y = gsp_filter(G,@(x) 1./(1+tau*x),randn(N,4));

%     fg.grad = @(x) 2*Mop(Mop(x)-y);
%     fg.eval = @(x) norm(Mop(x)-y)^2;
%     fg.beta = 2;
% %     paramtik.verbose = param.verbose -1;
% %     ftik.prox = @(x,T) gsp_prox_tik(x,tau * T,G,paramtik);
%     ftik.eval = @(x) tau* sum(gsp_norm_tik(G,x));
%     ftik.grad = @(x) 2 * G.L * x;
%     ftik.beta = 2 * G.lmax;
    param.maxit = 1000;
    param.tol = 1e-13;
    sol1 = gsp_regression_tik(G ,M, y , tau, param );
%    sol2 = solvep(y,{fg,ftik}, param );
    sol2 = (diag(M)+tau*G.L)^-1 * diag(M)*y;
    
    if norm(sol1-sol2,'fro')/norm(sol1,'fro') < 1e-5
        fprintf('GSP_PROX_TIK: exact 2 ok\n');
    else
        warning('GSP_PROX_TIK: exact 2 Pas OK');
        badness = norm(sol1-sol2,'fro')/norm(sol1,'fro') 
        errors = errors+1;
    end
%%

end



function errors = test_knn_classify_graph()

errors = 0;

% [ x,y,xx,yy,f,ff ] = prepare_usps_full( );
[x, y, xx, yy] = load_usps_full();

param.k = 6;
param.use_flann = 1;
gsp_knn_classify_graph( x(:,3:end)', xx', param );


param.k = 6;
param.use_flann = 0;
G = gsp_knn_classify_graph( x(:,3:end)', xx', param );

param.weighted = 1;
G2 = gsp_knn_classify_graph( x(:,3:end)', xx', param );


% Classification test

xtot = [x(:,3:end), xx];


M = zeros(size(xtot,2),1);
M(1:size(x(:,3:end),2)) = 1;
ytot = [y(3:end); yy];

paramr.tol = 0;
paramr.maxit = 1000;

sol1 = gsp_regression_knn(G,M,ytot);
sol2 = gsp_regression_tik(G,M,ytot,0,paramr);

if norm(sol1-sol2)/norm(sol1) < 1e-10
   fprintf('TEST_GRAPH_ML: gsp_kkn_classify_graph ok\n');
else
    errors = errors +1;
    norm(sol1-sol2)/norm(sol1)
    warning('TEST_GRAPH_ML: gsp_kkn_classify_graph pas OK!')
end


sol1 = gsp_classification_knn(G,M,ytot);
sol2 = gsp_classification_tik(G,M,ytot,0,paramr);

if norm(sol1-sol2)/norm(sol1) < 0.05
   fprintf('TEST_GRAPH_ML: gsp_kkn_classify_graph 2 ok\n');
else
    errors = errors +1;
    norm(sol1-sol2)/norm(sol1)
    warning('TEST_GRAPH_ML: gsp_kkn_classify_graph 2 pas OK!')
end





sol1 = gsp_regression_knn(G2,M,ytot);

sol2 = gsp_regression_tik(G2,M,ytot,0,paramr);

if norm(sol1-sol2)/norm(sol1) < 1e-4
   fprintf('TEST_GRAPH_ML: gsp_kkn_classify_graph 3 ok\n');
else
    errors = errors +1;
    norm(sol1-sol2)/norm(sol1)
    warning('TEST_GRAPH_ML: gsp_kkn_classify_graph 3 pas OK!')
end


sol1 = gsp_classification_knn(G2,M,ytot);

sol2 = gsp_classification_tik(G2,M,ytot,0,paramr);

if norm(sol1-sol2)/norm(sol1) < 1e-6
   fprintf('TEST_GRAPH_ML: gsp_kkn_classify_graph 4 ok\n');
else
    errors = errors +1;
    norm(sol1-sol2)/norm(sol1)
    warning('TEST_GRAPH_ML: gsp_kkn_classify_graph 4 pas OK!')
end



end


function errors = compare_flann_matlab()

    errors = 0;

    [x, y, xx, yy] = load_usps_full();


    param.k = 6;
    param.use_flann = 1;
    G1 = gsp_knn_classify_graph( x(:,3:end)', xx', param );


    param.use_flann = 0;
    G2 = gsp_knn_classify_graph( x(:,3:end)', xx', param );
    
    E = full(sum(sum((G1.W -G2.W)>0))/2);
    
    fprintf('Number of different edges: %d ', E)
    errors = errors + E>0;
    warning('This test is strange!!!!')
    
end






