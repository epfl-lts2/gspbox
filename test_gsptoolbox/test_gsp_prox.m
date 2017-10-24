function [ errors ] = test_gsp_prox( )

errors = 0;

errors = errors + test_norm_tik();

errors = errors + test_prox_tv();
errors = errors + test_norm_tv();
errors = errors + test_prox_tv_mat();
% errors = errors + test_prox_tv2(); % to be done
errors = errors + test_norm_tv2();
errors = errors + test_prox_tik();
errors = errors + test_prox_tik_mat();
errors = errors + test_prox_tik_A();



errors = errors + test_prox_tvmultidim();
errors = errors + test_prox_tikmultidim();



end


function errors = test_norm_tik()
errors = 0;

Ng = 10;
k = 4;
x = 5*rand(Ng,k);


G = gsp_sensor(Ng);

y1 = trace(x'*G.L*x);
y2 =sum(gsp_norm_tik(G,x));


if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-10
   errors = errors +1;
   warning('OPERATOR: Error in the norm tik test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: norm tik 1 ok\n');
end



y1 =gsp_norm_tik(G,x);
y2 = zeros(k,1);
for ii = 1:k
    y2(ii)=gsp_norm_tik(G,x(:,ii));
end
    
if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-10
   errors = errors +1;
   warning('OPERATOR: Error in the norm tik test 2')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: norm tik 2ok\n');
end




end



function errors = test_norm_tv()
errors = 0;

Ng = 10;
x = 5*rand(Ng,1);


G = gsp_path(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

y1 = norm_tv1d(x);
y2 = gsp_norm_tv(G,x(:));


if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-5
   errors = errors +1;
   warning('OPERATOR: Error in the norm_tv test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: norm_tv 1 ok\n');
end

k = 4;

x = 5*rand(Ng,k);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);


y1 =gsp_norm_tv(G,x);
y2 = zeros(k,1);
for ii = 1:k
    y2(ii)=gsp_norm_tv(G,x(:,ii));
end
    
if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-10
   errors = errors +1;
   warning('OPERATOR: Error in the norm tv test 2')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: norm tv 2ok\n');
end



end


function errors = test_prox_tv()
errors = 0;

Ng = 10;
x = 5*rand(Ng,1);


G = gsp_path(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param_tv.verbose = 0;
param_tv.tol = eps;
param_tv.maxit = 500;
y1 = gsp_prox_tv(x,1,G,param_tv);

y2 = prox_tv1d(x,1,param_tv);


if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-5
   errors = errors +1;
   warning('OPERATOR: Error in the prox_tv test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox_tv 1 ok\n');
end



end


function errors = test_prox_tv_mat()
errors = 0;

Ng = 500;
x = 5*rand(Ng,1);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param_tv.verbose = 0;
param_tv.tol = eps;
param_tv.maxit = 500;
param_tv.use_matrix = 0;
tic;
y1 = gsp_prox_tv(x,1,G,param_tv);
time_op = toc;
param_tv.use_matrix = 1;
tic;
y2 = gsp_prox_tv(x,1,G,param_tv);
time_mat = toc;

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-5
   errors = errors +1;
   warning('OPERATOR: Error in the prox_tv_mat test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox_tv mat 1 ok\n');
end



end


function errors = test_norm_tv2()
errors = 0;

Ng = 10;
x = 5*rand(Ng);


G = gsp_2dgrid(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

y1 = sum(norm_tv1d(x))+sum(norm_tv1d(x'));
y2 = gsp_norm_tv(G,x(:));


if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-5
   errors = errors +1;
   warning('OPERATOR: Error in the norm_tv test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: norm_tv 1 ok\n');
end



end


% function errors = test_prox_tv2()
% errors = 0;
% 
% Ng = 10;
% x = 5*rand(Ng);
% 
% 
% G = gsp_2dgrid(Ng);
% 
% G = gsp_estimate_lmax(G);
% G = gsp_adj2vec(G);
% 
% param_tv.verbose = 0;
% param_tv.tol = eps;
% param_tv.maxit = 200;
% y1 = gsp_prox_tv(x(:),1,G,param_tv);
% 
% y2 = prox_tv(x,1,param_tv);
% 
% 
% if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-5
%    errors = errors +1;
%    warning('OPERATOR: Error in the prox_tv test 1')
%    norm(y1(:)-y2(:))/norm(y1(:))
% 
% else
%     fprintf('OPERATOR: prox_tv 1 ok\n');
% end
% 
% 
% 
% end


function errors = test_prox_tik()
errors = 0;

Ng = 500;
x = 5*rand(Ng,1);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param.verbose = 0;
tic;
param.order = 100;
y1 = gsp_prox_tik(x,1,G,param);
time_cheb = toc;
param.A = @(x) x;
param.At = @(x) x;
param.tol = eps(1000);
param.maxit = 500;
tic;
param.use_matrix = 1;
y2 = gsp_prox_tik(x,1,G,param);
time_pcg = toc;
if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-6
   errors = errors +1;
   warning('OPERATOR: Error in the prox tik test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tik 1 ok\n');
end



end


function errors = test_prox_tik_mat()
errors = 0;

Ng = 500;
x = 5*rand(Ng,1);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param.verbose = 0;


param.A = @(x) x;
param.At = @(x) x;
param.tol = eps(1000);
param.maxit = 500;

param.use_matrix = 0;
tic;
y1 = gsp_prox_tik(x,1,G,param);
time_grad = toc;
param.use_matrix = 1;
tic;
y2 = gsp_prox_tik(x,1,G,param);
time_mat = toc;

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-4
   errors = errors +1;
   warning('OPERATOR: Error in the prox tik mat test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tik mat 1 ok\n');
end



end



function errors = test_prox_tvmultidim()
errors = 0;

Ng = 30;
k = 10;
x = 5*rand(Ng,k);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param.verbose = 0;
param.maxit = 100;
param.tol = eps;
y1 = gsp_prox_tv(x,1,G,param);

y2 = zeros(Ng,k);
for ii = 1:k
    y2(:,ii) = gsp_prox_tv(x(:,ii),1,G,param);
end

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-6
   errors = errors +1;
   warning('OPERATOR: Error in the prox tv multi dim test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tv multi dim ok\n');
end


end



function errors = test_prox_tikmultidim()
errors = 0;

Ng = 30;
k = 10;
x = 5*rand(Ng,k);


G = gsp_sensor(Ng);

G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

param.verbose = 0;
param.maxit = 100;
param.tol = 1000*eps;
y1 = gsp_prox_tik(x,1,G,param);

y2 = zeros(Ng,k);
for ii = 1:k
    y2(:,ii) = gsp_prox_tik(x(:,ii),1,G,param);
end

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-6
   errors = errors +1;
   warning('OPERATOR: Error in the prox tik multi dim test 1')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tik multi dim 1 ok\n');
end


param.A = @(x) x;

param.At = @(x) x;

y1 = gsp_prox_tik(x,1,G,param);

y2 = zeros(Ng,k);
for ii = 1:k
    y2(:,ii) = gsp_prox_tik(x(:,ii),1,G,param);
end

if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-6
   errors = errors +1;
   warning('OPERATOR: Error in the prox tik multi dim test 2')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tik multi dim 2 ok\n');
end


end



function errors = test_prox_tik_A()
errors = 0;

Ng = 50;
x = 5*rand(Ng,1);

A = rand(Ng);

G = gsp_sensor(Ng);

G2 = G;


G = gsp_estimate_lmax(G);
G = gsp_adj2vec(G);

Lp = A'*G.L'*A;
G2.L = Lp;
G2 = gsp_estimate_lmax(G2);
param.verbose = 0;
param.order = 100;
y1 = gsp_prox_tik(x,1,G2,param);
param.A = @(x) A*x;
param.At = @(x) A'*x;
param.tol = eps(10000);
param.maxit = 500;
param.use_matrix = 1;
y2 = gsp_prox_tik(x,1,G,param);
if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-4
   errors = errors +1;
   warning('OPERATOR: Error in the prox tik A')
   norm(y1(:)-y2(:))/norm(y1(:))

else
    fprintf('OPERATOR: prox tik A ok\n');
end



end
