function [ errors ] = test_operators( )

errors = 0;

errors = errors + test_fourier_basis_ring();
errors = errors + test_fourier_basis_ring_advanced();

%errors = errors + test_translate();
errors = errors + test_modulate();

errors = errors + test_gwft();
errors = errors + test_gft();
errors = errors + test_ngwft();

errors = errors + test_divgradient();
errors = errors + test_gradient();
errors = errors + test_gradient2();
errors = errors + test_gradient_mat();

errors = errors + test_directed();
errors = errors + test_advection();

errors = errors + test_div();
errors = errors + test_div2();
errors = errors + test_div3();
errors = errors + test_laplacian();

errors = errors + test_kron_reduction();



end

function s = eigvec_ring(input,N)
s = zeros(N,numel(input));
n = (0:(N-1))';
for ind = 1:numel(input)
    if mod(input(ind),2) % odd
        s(:,ind) = 1/sqrt(N)*exp(2*pi*1i*(N-(input(ind)+1)/2)*n/N);
    else % even     
        s(:,ind) = 1/sqrt(N)*exp(pi*1i*input(ind)*n/N);    
    end
end

end

function s = eigval_ring(input,N)
s = zeros(size(input));
for ind = 1:numel(input)
    if mod(input(ind),2) % odd
        s(ind) = 1-cos(pi * (input(ind)+1)/N );
    else % even     
        s(ind) = 1-cos(pi * (input(ind))/N );
    end
end

end

function errors = test_fourier_basis_ring_advanced()
errors = 0;
N = 10;
e = 2*eigval_ring(0:N-1,N)';
U = eigvec_ring(0:N-1,N);
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);
errors = errors +gsp_assert_test(G.e,e,1e-14,'test ring eigenvalues N=10');
errors = errors +gsp_assert_test(G.U,U,1e-14,'test ring eigenvector N=10');

N = 11;
e = 2*eigval_ring(0:N-1,N)';
U = eigvec_ring(0:N-1,N);
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);
errors= errors +gsp_assert_test(G.e,e,1e-14,'test ring eigenvalues N=11');
errors = errors +gsp_assert_test(G.U,U,1e-14,'test ring eigenvector N=11');

end


function errors = test_laplacian()
errors = 0;
N = 10;
G = gsp_path(N);


x = rand(N,1);

y1 = G.L*x;
y2 = -laplacianx_op(x);



if norm(y1-y2)/norm(y1)>10e-10;
   errors = errors +1;
   norm(y1-y2)/norm(y1)
   warning('OPERATOR: Error in laplacian test 1')

else
    fprintf('OPERATOR: laplacian 1 ok\n');
end


end


function errors = test_div()
errors = 0;
N = 10;
G = gsp_path(N);
G = gsp_adj2vec(G);

x = rand(N-1,1);

y1 = gsp_div(G,x);
y2 = -div_op1d([x;0]);



if norm(y1-y2)/norm(y1)>10e-10;
   errors = errors +1;
   norm(y1-y2)/norm(y1)
   warning('OPERATOR: Error in div test 1')

else
    fprintf('OPERATOR: div 1 ok\n');
end


end



function errors = test_div2()
errors = 0;
N = 10;
G = gsp_sensor(N);
G = gsp_adj2vec(G);

M = 5;

x = rand(N,M);

t = gsp_grad(G,x);

y1 = gsp_div(G,t);
y2 = gsp_div_old(G,t);



if norm(y1(:)-y2(:))/norm(y1(:))>10e-10;
   errors = errors +1;
   norm(y1(:)-y2(:))/norm(y1(:))
   warning('OPERATOR: Error in div test 2')

else
    fprintf('OPERATOR: div 2 ok\n');
end


end

function errors = test_div3()
errors = 0;
N = 10;
G = gsp_sensor(N);
G = gsp_adj2vec(G);

M = 5;

x = rand(N,M);

t = gsp_grad(G,x);

y1 = gsp_div(G,t);

y2 = zeros(N,M);
for ii = 1:M
    y2(:,ii) = gsp_div_old(G,t(:,ii));
end



if norm(y1(:)-y2(:))/norm(y1(:))>10e-10;
   errors = errors +1;
   norm(y1(:)-y2(:))/norm(y1(:))
   warning('OPERATOR: Error in div test 3')

else
    fprintf('OPERATOR: div 3 ok\n');
end


end



function errors = test_gradient()
errors = 0;
N = 50;
G = gsp_path(N);
G = gsp_adj2vec(G);

x = rand(N,1);

y1 = gsp_grad(G,x);
y2 = gradient_op1d(x);



if norm(y1-y2(1:(end-1)))/norm(y1)>10e-10;
   errors = errors +1;
   norm(y1-y2(1:(end-1)))/norm(y1)
   warning('OPERATOR: Error in gradient test 1')

else
    fprintf('OPERATOR: gradient 1 ok\n');
end


end

function errors = test_gradient2()
errors = 0;
N = 50;
G = gsp_sensor(N);
G = gsp_adj2vec(G);

M = 5;

x = rand(N,M);

y1 = gsp_grad(G,x);

y2 = zeros(length(G.v_in),M);
for ii = 1:M
    y2(:,ii) = gsp_grad(G,x(:,ii));
end



if norm(y1(:)-y2(:))/norm(y1(:))>10e-10;
   errors = errors +1;
   norm(y1(:)-y2(:))/norm(y1(:))
   warning('OPERATOR: Error in gradient test 2')

else
    fprintf('OPERATOR: gradient 2 ok\n');
end


end

function errors = test_gradient_mat()
errors = 0;
N = 50;
G = gsp_sensor(N);
G = gsp_adj2vec(G);


x = rand(N,1);

y1 = gsp_grad(G,x);

D = gsp_grad_mat(G);


y2 = D*x;



if norm(y1(:)-y2(:))/norm(y1(:))>10e-10;
   errors = errors +1;
   norm(y1(:)-y2(:))/norm(y1(:))
   warning('OPERATOR: Error in gradient mat test')

else
    fprintf('OPERATOR: gradient mat ok\n');
end


end


function errors = test_directed()
errors = 0;

try
    
    G = gsp_sensor(50);
    
    Gd = gsp_assign_rand_direction(G,.5);
    Gd = gsp_adj2vec(Gd);
    
    fprintf('OPERATOR: directed ok\n');
catch
    errors = errors + 1;
    warning('OPERATOR: Error directed test')
end

end

function errors = test_advection()
errors = 0;
G = gsp_sensor();

Gd = gsp_assign_rand_direction(G,.5);
Gd = gsp_adj2vec(Gd);


W = diag(diag(Gd.Adv)) - Gd.Adv';


err = norm(Gd.W-W,'fro')/norm(Gd.W,'fro');

if err>1e-10;
   errors = errors +1;
   err
   warning('OPERATOR: Error advection test')

else
    fprintf('OPERATOR: advection ok\n');
end

end





function errors = test_divgradient()
errors = 0;
N = 50;
G = gsp_sensor(N);
G = gsp_adj2vec(G);

x = rand(N,1);

y1 = gsp_div(G,gsp_grad(G,x));
y2 = G.L*x;



if norm(y1-y2)/norm(y1)>10e-10
   errors = errors +1;
   norm(y1-y2)/norm(y1)
   warning('OPERATOR: Error in divgradient test 1')

else
    fprintf('OPERATOR: divgradient 1 ok\n');
end

if norm(full(G.Diff'*G.Diff-G.L))>10e-10;
   errors = errors +1;
   norm(full(G.Diff'*G.Diff-G.L))
   warning('OPERATOR: Error in divgradient test 2')

else
    fprintf('OPERATOR: divgradient 2 ok\n');
end


G = gsp_create_laplacian(G,'normalized');
G = gsp_adj2vec(G);

y1 = gsp_div(G,gsp_grad(G,x));
y2 = G.L*x;



if norm(y1-y2)/norm(y1)>10e-10
   errors = errors +1;
   norm(y1-y2)/norm(y1)
   warning('OPERATOR: Error in divgradient test 3')

else
    fprintf('OPERATOR: divgradient 3 ok\n');
end

if norm(full(G.Diff'*G.Diff-G.L))>10e-10;
   errors = errors +1;
   norm(full(G.Diff'*G.Diff-G.L))
   warning('OPERATOR: Error in divgradient test 3')

else
    fprintf('OPERATOR: divgradient 3 ok\n');
end


end

function errors = test_fourier_basis_ring()
errors = 0;
N = 50;
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);



if sum(diff(G.e) < -1e-10)
   errors = errors +1;
   warning('OPERATOR: Error in the Fourier test 1')

else
    fprintf('OPERATOR: Fourier 1 ok\n');
end

if sum(norm(G.L-G.U*diag(G.e)*G.U','fro') < -1e-10)
   errors = errors +1;
   warning('OPERATOR: Error in the Fourier test 2')

else
    fprintf('OPERATOR: Fourier 2 ok\n');
end

N = 51;
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);

if sum(norm(G.L-G.U*diag(G.e)*G.U','fro') < -1e-10)
   errors = errors +1;
   warning('OPERATOR: Error in the Fourier test 2')

else
    fprintf('OPERATOR: Fourier 2 ok\n');
end


end

function errors = test_translate()
errors = 0;
N = 50;
k = 10;
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);

f = rand(N,1);

ft = gsp_translate(G,f,k);

if norm(ft-circshift(f,[-k+1,0]))<1e-10
    fprintf('OPERATOR: translate ok\n');
else
    errors = errors +1;
   warning('OPERATOR: Error in the translate test')

end

end

function errors = test_modulate()
errors = 0;
N = 30;
k = 10;
G = gsp_ring(N);
G = gsp_compute_fourier_basis(G);

f = rand(N,1);

% compute the ordering vectors
v = gsp_classic2graph_eig_order( N );
iv = 1:N;
iv(v) = iv;

ft = gsp_modulate(G,f,iv(k));

fthat(v) = gsp_gft(G,ft);
ft2hat(v) = gsp_gft(G,f); % back into the classical ordering
ft2hat = circshift(ft2hat,[0,-k+1]); % circshift


if norm(fthat-ft2hat)<1e-10
    fprintf('OPERATOR: modulate ok\n');
else
    errors = errors +1;
   warning('OPERATOR: Error in the modulate test')

end

end


function errors = test_kron_reduction()

errors = 0;

try
    figure(100)
      N = 64;
      param.distribute = 1;
      param.Nc = 5;
      param.regular = 1;
      G = gsp_sensor(N,param);
      ind = 1:2:N;
      Gnew = gsp_kron_reduce( G,ind );
      subplot(121)
      gsp_plot_graph(G);
      title('Original graph');
      subplot(122)
      gsp_plot_graph(Gnew);
      title('Kron reduction');

    close(100)

    fprintf('OPERATOR: KRON ok\n');
catch
    errors = errors +1;
   warning('OPERATOR: Error in the KRON test')

end

gsp_reset_seed(1)
      G1 = gsp_kron_reduce( G,ind );
      gsp_reset_seed(1)
      Lnew= gsp_kron_reduce_old( G.L,ind );
if norm(Lnew-G1.L,'fro')<1e-10
      fprintf('OPERATOR: KRON ok\n');
else
          errors = errors +1;
    norm(G2.W-G1.W,'fro')
   warning('OPERATOR: Error in the KRON test')
end

end


function errors = test_gft()

N = 30;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

X = rand(N,20,10);

Xhat1 = gsp_gft(G,X);

Xhat2 = zeros(N,20,10);
for ii = 1:10
   Xhat2(:,:,ii) = gsp_gft(G, X(:,:,ii)); 
end

errors = gsp_assert_test(Xhat1,Xhat2,1e-10,'GFT');

end