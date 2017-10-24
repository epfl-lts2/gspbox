function errors = test_gsp_distanz()
    

errors = 0;
errors = errors + test_flann();

errors = errors + test_full();

errors = errors + test_use_full();

end

function errors = test_use_full()

errors = 0;

param.k = 4;
param.epsilon = 4;
n1 = 20;
n2 = 60;
X1 = rand(100,n1);
X2 = rand(100,n2);

param.use_full = 0;
param.type = 'knn';

[indx, indy, dist] = gsp_nn_distanz(X1,X2,param);
A1 = sparse(indx,indy,dist,n1,n2);

param.use_full = 1;
[indx, indy, dist] = gsp_nn_distanz(X1,X2,param);
A2 = sparse(indx,indy,dist,n1,n2);


if norm(A1 - A2,'fro')/norm(A1 ,'fro') > 1e-10
   norm(A1 - A2,'fro')/norm(A1 ,'fro')
   errors = errors +1;
   warning('DISTANZ: Error in test use_full knn')
figure    
subplot(121)
imagesc(A1)
subplot(122)
imagesc(A2)
else
   fprintf('DISTANZ: test use_full knn ok\n');

   
end

param.use_full = 0;
param.type = 'radius';

[indx, indy, dist] = gsp_nn_distanz(X1,X2,param);
A1 = sparse(indx,indy,dist,n1,n2);

param.use_full = 1;
[indx, indy, dist] = gsp_nn_distanz(X1,X2,param);
A2 = sparse(indx,indy,dist,n1,n2);


if norm(A1 - A2,'fro')/norm(A1 ,'fro') > 1e-10
   norm(A1 - A2,'fro')/norm(A1 ,'fro')
   errors = errors +1;
   warning('DISTANZ: Error in test use_full radius')
figure    
subplot(121)
imagesc(A1)
subplot(122)
imagesc(A2)
else
   fprintf('DISTANZ: test use_full radius ok\n');

end


end


function errors = test_flann()
%gsp_reset_seed

%test flan l_2 
errors = 0;
N = 10;
Nfeat = 5;
x1 = rand(Nfeat,N);
x2 = rand(Nfeat,N);

param.k = N;
param.rescale = 0;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x1,x2,param); 
y1 = full(sparse(indx,indy,dist,N,N));

param.use_flann = 1;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x1,x2,param); 

y2 = full(sparse(indx,indy,dist,N,N));

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test flann 1')
   norm(y1(:)-y2(:))/norm(y1(:))
   figure
   subplot(121)
   imagesc(y1)
   caxis([min(y1(:)),max(y1(:))])
   subplot(122)
   imagesc(y2)
  caxis([min(y1(:)),max(y1(:))])


else
    fprintf('DISTANZ: test flann 1 ok\n');
end

%test flan l_1 
param.use_l1 = 1;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x1,x2,param); 
y1 = full(sparse(indx,indy,dist,N,N));

param.use_flann = 1;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x1,x2,param); 
y2 = full(sparse(indx,indy,dist,N,N));



if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test flann 2')
   norm(y1(:)-y2(:))/norm(y1(:))
   figure
   subplot(121)
   imagesc(y1)
   caxis([min(y1(:)),max(y1(:))])
   subplot(122)
   imagesc(y2)
  caxis([min(y1(:)),max(y1(:))])


else
    fprintf('DISTANZ: test flann 2 ok\n');
end



%%
%test flann l_2 
param.use_l1 = 0;
param.use_flann = 0;
param.sigma = 0.1;
x = rand(1000,10);
G0 = gsp_nn_graph(x,param);
param.use_flann = 1;
param.flann_checks = 10000000;
t1 = tic ;
G1 = gsp_nn_graph(x,param);
time1 = toc(t1);
param.use_flann = 1;
param.flann_checks = 256;
t2 = tic;
G2 = gsp_nn_graph(x,param);
time2 = toc(t2);

if norm(G1.W-G0.W,'fro')/norm(G0.W,'fro') > 0.01
   norm(G1.W-G0.W,'fro')/norm(G0.W,'fro')
   errors = errors +1;
   warning('DISTANZ: Error in test flann 3')
else
   fprintf('DISTANZ: test flann 3 ok\n');

end


%test flann l_1 
param.use_l1 = 1;
param.use_flann = 0;
param.sigma = 0.1;
x = rand(1000,10);
G0 = gsp_nn_graph(x,param);
param.use_flann = 1;
param.flann_checks = 10000000;
t1 = tic ;
G1 = gsp_nn_graph(x,param);
time1 = toc(t1);
param.use_flann = 1;
param.flann_checks = 256;
t2 = tic;
G2 = gsp_nn_graph(x,param);
time2 = toc(t2);

if norm(G1.W-G0.W,'fro')/norm(G0.W,'fro') > 0.01
   norm(G1.W-G0.W,'fro')/norm(G0.W,'fro')
   errors = errors +1;
   warning('DISTANZ: Error in test flann 4')
else
   fprintf('DISTANZ: test flann 4 ok\n');

end
% norm(G1.W-G0.W,'fro')/norm(G0.W,'fro')
% norm(G2.W-G0.W,'fro')/norm(G0.W,'fro')


end


function errors = test_full()

errors = 0;
N = 10;
N2 = 20;
Nfeat = 5;
x = rand(Nfeat,N);

y1 = gsp_distanz(x);
param.k = N;
param.rescale = 0;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x,x,param); 

y2 = full(sparse(indx,indy,dist,N,N));

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test full 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('DISTANZ: test full 1 ok\n');
end

x2 = rand(Nfeat,N2);
y1 = gsp_distanz(x,x2);
param.k = N;
param.rescale = 0;
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x,x2,param); 
y2 = full(sparse(indx,indy,dist,N,N2));



[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x,x2,param); 
y2 = full(sparse(indx,indy,dist,N,N2));

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test full 2')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('DISTANZ: test full 2 ok\n');
end

y1 = gsp_distanz(x);
param.type = 'radius';
param.rescale = 0;
param.epsilon = max(y1(:)+0.1);
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x,x,param); 
y2 = full(sparse(indx,indy,dist,N,N));

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test full 3')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('DISTANZ: test full 3 ok\n');
end

y1 = gsp_distanz(x,x2);
param.type = 'radius';
param.rescale = 0;
param.epsilon = max(y1(:)+0.1);
[indx, indy, dist, Xo1, Xo2, epsilon] = gsp_nn_distanz(x,x2,param); 
y2 = full(sparse(indx,indy,dist,N,N2));
if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-8
   errors = errors +1;
   warning('DISTANZ: Error in test full 4')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('DISTANZ: test full 4 ok\n');
end

end




