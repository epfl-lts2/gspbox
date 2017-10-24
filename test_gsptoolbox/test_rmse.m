function errors = test_rmse()
errors = 0;

errors = errors + test_rmse_mv();
errors = errors + test_nn_graph();


errors = errors + test_rmse_mv_graph();




end


function errors = test_rmse_mv()
errors = 0;

X = rand(10,100);

d1 = gsp_distanz(X);
d2 = gsp_rmse_mv(X)*sqrt(10);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST RMSE MV OK\n');
    
else
    fprintf('ERROR IN TEST RMSE MV !!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
end

end



function errors = test_nn_graph()
errors = 0;

X = rand(10,100);

sigma = 2;

d1 = exp(-gsp_distanz(X).^2/sigma);
d1 = d1-diag(diag(d1));
param.rescale = 0;
param.center = 0;
param.sigma = sigma;
param.type = 'radius';
param.epsilon = 10000;
G = gsp_nn_graph(X',param);
d2 = full(G.W);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST NN graph OK\n');
    
else
    fprintf('ERROR IN TEST NN graph !!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
end

X1 = X-repmat(mean(X'),100,1)';

d1 = exp(-gsp_distanz(X1).^2/sigma);
d1 = d1-diag(diag(d1));
param.rescale = 0;
param.center = 1;
param.sigma = sigma;
param.type = 'radius';
param.epsilon = 10000;
G = gsp_nn_graph(X',param);
d2 = full(G.W);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST NN graph 2 OK\n');
    
else
    fprintf('ERROR IN TEST NN graph 2 !!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
end

end





function errors = test_rmse_mv_graph()
errors = 0;

Nfeat = 3;
X = rand(Nfeat,100);

sigma = 2;


param.sigma = sigma;
param.type = 'radius';
param.epsilon = 10000;
param.rescale = 0;
param.center = 0;
W = gsp_rmse_mv(X)*sqrt(Nfeat);
W = exp(-W.^2/sigma);
W = W-diag(diag(W));

d1 = full(W);
G = gsp_rmse_mv_graph(X',param);
d2 = full(G.W);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST gsp_rmse_mv_graph 1 OK\n');
    
else
    fprintf('ERROR gsp_rmse_mv_graph  1!!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
end

param.rescale = 0;
param.epsilon = 0.5;
param.type = 'radius';
G1 = gsp_nn_graph(X',param);
G2 = gsp_rmse_mv_graph(X',param);
d1 = full(G1.W);
d2 = full(G2.W);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST gsp_rmse_mv_graph 2 OK\n');
    
else
    fprintf('ERROR gsp_rmse_mv_graph  2!!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
            figure()
    subplot(121)
    imagesc(G1.W)    
    subplot(122)
    imagesc(G2.W)

end

param.rescale = 0;
param.type = 'knn';
param.k = 5;
G1 = gsp_nn_graph(X',param);
G2 = gsp_rmse_mv_graph(X',param);
d1 = full(G1.W);
d2 = full(G2.W);

if norm(d1-d2,'fro')/norm(d1)<1e-8
    fprintf('TEST gsp_rmse_mv_graph 3 OK\n');
    
else
    fprintf('ERROR gsp_rmse_mv_graph  3!!!!!!!\n')
    norm(d1-d2,'fro')/norm(d1)  
    errors = errors +1;
        figure()
    subplot(121)
    imagesc(G1.W)    
    subplot(122)
    imagesc(G2.W)
end


end