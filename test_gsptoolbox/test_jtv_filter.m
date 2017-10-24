function errors = test_jtv_filter()
%TEST_JTV_FILTER This function test all the time-vertex filters
close all
errors = 0;
warning on

errors = errors + test_jtv_graph();
errors = errors + test_jtv_graph_diagonalization();

errors = errors + test_diffusion();
errors = errors + test_wave();
errors = errors + test_jft_ijft();
errors = errors + test_swap_transform();

errors = errors + test_theory()
errors = errors + test_swap_localization();


errors = errors + test_jtv_cheb_coeff();
errors = errors + test_jtv_cheb_op();

errors = errors + test_jtv_filter_evaluate();
errors = errors + test_jtv_filter_array();

errors = errors + test_jtv_filter_analysis();
errors = errors + test_jtv_filter_synthesis();

errors = errors + test_jtv_matrix_analysis_op();
errors = errors + test_jtv_matrix_synthesis_op();

errors = errors + test_jtv_filter_dual();
errors = errors + test_jtv_filter_inverse();



try  %#ok<TRYNC>
    close(100)
end

end

function errors = test_jtv_graph()
errors=0;

try
    
    N = 100;
    T = 200;
    fs = 1;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T,fs);
    
    fprintf('Test JTV - timevertex graph:  OK\n');
catch
    errors = errors + 1;
    warning('Test JTV - Error in timevertex graph test')
end

end

function errors = test_jtv_graph_diagonalization()
errors  = 0;
T=30;
N=20;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
DFT = dftmtx(T)';

Gt = gsp_ring(T);
Gp = gsp_graph_product(G,Gt);
Gp = gsp_create_laplacian(Gp,'combinatorial');
Gp = gsp_compute_fourier_basis(Gp);

U = kron(G.U,DFT);

errors = errors + gsp_assert_test(U'*Gp.L*U,diag(diag(U'*Gp.L*U)),1e-10,'JTV - diagonalization test');

end

function errors = test_diffusion()

errors = 0;
try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    [g,ft] = gsp_jtv_design_diffusion(G);
    gsp_plot_jtv_filter(G,g,ft);
    close
    
    fprintf('JTV: diffusion kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error diffusion kernel 1 test')
end

try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    tau=.1;
    [g,ft] = gsp_jtv_design_diffusion(G,tau);
    gsp_plot_jtv_filter(G,g,ft);
    close
    
    fprintf('JTV: diffusion kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error diffusion kernel 2 test')
end


try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    
    tau=[.01 0.1 1 2];
    param.normalize=1;
    param.show_sum=0;
    param.title=num2cell(tau);
    [g,ft] = gsp_jtv_design_diffusion(G,tau,param);
    
    gsp_plot_jtv_filter(G,g,ft);
    close
    
    
    fprintf('JTV: diffusion kernel 3 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error diffusion kernel 3 test')
end

end

function errors = test_wave()

errors = 0;
try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    
    [g,ft] = gsp_jtv_design_wave(G);
    gsp_plot_jtv_filter(G,g,ft);
    
    close
    
    fprintf('JTV: wave kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error wave kernel 1 test')
end

try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    
    alpha=1;
    [g,ft] = gsp_jtv_design_wave(G,alpha);
    gsp_plot_jtv_filter(G,g,ft);
    
    close
    
    
    fprintf('JTV: wave kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error wave kernel 2 test')
end


try
    
    N = 100;
    T = 200;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    
    alpha=1;
    param.normalize=1;
    [g,ft] = gsp_jtv_design_wave(G,alpha,param);
    gsp_plot_jtv_filter(G,g,ft);
    
    close
    
    
    fprintf('JTV: wave kernel 3 ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error wave kernel 3 test')
end

end


function errors = test_jft_ijft()
errors = 0;

try
    
    N=50;
    T=100;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_compute_fourier_basis(G);
    
    x = randn(N,T);
    
    xhat = gsp_jft(G,x);
    
catch
    errors = errors + 1;
    warning('JTV: Error jtgft test')
end

try
    
    
    N=50;
    T=100;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_compute_fourier_basis(G);
    
    xhat = randn(N,T);
    x = gsp_ijft(G,x);
    
    
catch
    errors = errors + 1;
    warning('JTV: Error ijtgft test')
end


N=50;
T=100;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

x = randn(N,T);

x2 = gsp_ijft(G,gsp_jft(G,x));



errors = errors + gsp_assert_test(x,x2,eps(1000),'JTV - inverse jtv Fourier transform');

x = rand(N,T,1);

s1 = gsp_jft(G,x);
s2 = gsp_jft_simple(G,x);

errors = errors + gsp_assert_test(s1,s2,eps(1000),'JTV - JFT 1 dim');


x = rand(N,T,10);

s1 = gsp_jft(G,x);
s2 = gsp_jft_simple(G,x);

errors = errors + gsp_assert_test(s1,s2,eps(1000),'JTV - JFT 2 dim');


end

function errors = test_swap_transform()
errors = 0;

N=50;
T=100;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

x = randn(N,T);

x1 = fft(gsp_gft(G,x),[],2)/sqrt(T);
x2 = gsp_gft(G,fft(x,[],2))/sqrt(T);
x3 = gsp_jft(G,x);

e1 = norm(x1-x2,'fro')/norm(x1,'fro');
e2 = norm(x1-x3,'fro')/norm(x1,'fro');
e3 = norm(x2-x3,'fro')/norm(x2,'fro');

errors = errors + gsp_assert_test(0,mean([e1 e2 e3]),eps(1000),'JTV - commutative fft gft transform');

%

xhat = gsp_jft(G,x);

x1 = ifft(gsp_igft(G,xhat),[],2)*sqrt(T);
x2 = gsp_igft(G,ifft(xhat,[],2))*sqrt(T);
x3 = gsp_ijft(G,xhat);

e1 = norm(x1-x2,'fro')/norm(x1,'fro');
e2 = norm(x1-x3,'fro')/norm(x1,'fro');
e3 = norm(x2-x3,'fro')/norm(x2,'fro');

errors = errors + gsp_assert_test(0,mean([e1 e2 e3]),eps(1000),'JTV - commutative ifft igft transform');


end


function errors = test_swap_localization()
errors = 0;

N = 100;
T = 50;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

x = rand(N,T);
vertex=randi(N);
time=randi(T);

% to be updated (with gsp_localize?)
error('Does the boundary still make sense?')
param.boundary='periodic';
x1p = gsp_time_translate(G,gsp_translate_old(G,x,vertex),time,param);
x2p = gsp_translate_old(G,gsp_time_translate(G,x,time,param),vertex);
param.boundary='symmetric';
x1s = gsp_time_translate(G,gsp_translate_old(G,x,vertex),time,param);
x2s = gsp_translate_old(G,gsp_time_translate(G,x,time,param),vertex);
param.boundary='absorbing';
x1a = gsp_time_translate(G,gsp_translate_old(G,x,vertex),time,param);
x2a = gsp_translate_old(G,gsp_time_translate(G,x,time,param),vertex);

n1=norm(x1p-x2p)/norm(x1p);
n2=norm(x1s-x2s)/norm(x1s);
n3=norm(x1a-x2a)/norm(x1a);

if  (n1+n2+n3)/3< eps(1000)
    fprintf('JTV: commutative translation operator ok\n');
    
else
    errors = errors + 1;
    warning('JTV: error commutative translation operator test')
end






end

function errors = test_jtv_cheb_coeff()
errors = 0;

try
    N = 100;
    T = 50;
    alpha = 1;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    [g,ft] = gsp_jtv_design_wave(G, alpha);
    c = gsp_jtv_cheby_coeff(G, g,ft);
    
    fprintf('FILTER: jtv cheb coeff 1 ok\n');
    
catch
    errors = 1;
    warning('FILTER: Error in jtv cheb coeff 1 test')
end

try
    N = 100;
    T = 50;
    alpha = [ 0.1 0.5 1 ];
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    [g,ft] = gsp_jtv_design_wave(G, alpha);
    c = gsp_jtv_cheby_coeff(G, g,ft);
    
    fprintf('FILTER: jtv cheb coeff 2 ok\n');
    
catch
    errors = 1;
    warning('FILTER: Error in jtv cheb coeff 2 test')
end



end

function errors = test_jtv_cheb_op()

errors = 0;
try
    N = 100;
    T = 50;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_estimate_lmax(G);
    alpha = 1;
    [g,ft] = gsp_jtv_design_wave(G, alpha);
    c = gsp_jtv_cheby_coeff(G, g,ft);
    f = randn(N,T);
    r = gsp_jtv_cheby_op(G, c, f);
    
    fprintf('FILTER: cheb op ok\n');
catch
    errors = 1;
    warning('FILTER: Error in cheb op test')
end

end

function errors = test_jtv_filter_evaluate()
errors = 0;

N = 100;
T = 100;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
G = gsp_jtv_graph(G,T);

alpha=1;
[g,ft] = gsp_jtv_design_wave(G,alpha);
x = G.e;
param.filtertype = ft;
t = gsp_jtv_ta(G);
f = gsp_jtv_fa(G);

x1 = gsp_jtv_filter_evaluate_simple(G,g,x,param);

x2 = gsp_jtv_filter_evaluate(g,ft,x,t,param);

errors = errors + gsp_assert_test(x1,x2,1e-12,'JTV filter evaluate 1');

g = @(x,y) x+y;
ft = 'js';
param.filtertype = ft;
x1 = gsp_jtv_filter_evaluate_simple(G,g,x,param);

x2 = gsp_jtv_filter_evaluate(g,ft,x,f);

errors = errors + gsp_assert_test(x1,x2,1e-12,'JTV filter evaluate 2');

end



function errors = test_jtv_filter_array()
errors = 0;

N = 500;
T = 100;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
G = gsp_jtv_graph(G,T);

alpha=[0 0.5 0.7 1];
[g,ftg] = gsp_jtv_design_wave(G,alpha);
x = G.e;
t = gsp_jtv_ta(G);

[h,fth] = gsp_jtv_filter_array(G,g,ftg);

x1 = gsp_jtv_filter_evaluate(g,ftg,x,t);
x2 = gsp_jtv_filter_evaluate(h,fth,x,t);


errors = errors + gsp_assert_test(x1,x2,1e-12,'JTV filter array 1');


g = @(x,y) exp((x-y.^2)/0.1);
ftg = 'js';
f = gsp_jtv_fa(G);

[h,fth] = gsp_jtv_filter_array(G,g,ftg);

x1 = gsp_jtv_filter_evaluate(g,ftg,x,f);
x2 = gsp_jtv_filter_evaluate(h,fth,x,f);


errors = errors + gsp_assert_test(x1,x2,1e-12,'JTV filter array 2');

end



function errors = test_jtv_filter_analysis()
%%
errors = 0;

try
    N = 100;
    T = 500;
    G = gsp_sensor(N);
    G = gsp_jtv_graph(G,T);
    G = gsp_compute_fourier_basis(G);
    
    alpha=[0.1:0.1:0.7];
    [g,ft] = gsp_jtv_design_wave(G,alpha);
    
    x = rand(N,T);
    
    param.method = 'exact';
    coeff = gsp_jtv_filter_analysis(G,g,ft,x,param);
    
    fprintf('JTV: analysis ok\n');
catch
    errors = errors + 1;
    warning('JTV: Error analysis test')
end


N = 1000;
T = 100;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_estimate_lmax(G);

alpha=[0.1:0.1:0.7];
[g,ft] = gsp_jtv_design_wave(G,alpha);

x = randn(N,T);

param.method = 'exact';
tic
G = gsp_compute_fourier_basis(G);
c_exact = gsp_jtv_filter_analysis(G,g,ft,x,param);
param.filtertype = ft;
c_exact2 = gsp_jtv_filter_analysis_simple(G,g,x,param);
time=toc;
fprintf(['Eigendecomposition + Exact filtering time: ' num2str(time) '\n'])

errors = errors + gsp_assert_test(c_exact,c_exact2,1e-10,'JTV - exact analysis 1');


param.method = 'cheby';
param.order=40;
tic
c_cheby = gsp_jtv_filter_analysis(G,g,ft,x,param);
time=toc;
fprintf(['Cheby filtering time: ' num2str(time) '\n'])


errors = errors + gsp_assert_test(c_exact,c_cheby,1e-4,'JTV - cheby analysis 1');


param.method = 'exact';
g = @(x,y) exp(-x.*(y.^2)/0.1);
ft = 'js';
c_exact = gsp_jtv_filter_analysis(G,g,ft,x,param);

param.filtertype = ft;
c_exact2 = gsp_jtv_filter_analysis_simple(G,g,x,param);
errors = errors + gsp_assert_test(c_exact,c_exact2,1e-10,'JTV - exact analysis 2');


param.method = 'cheby';
param.order=20;
c_cheby = gsp_jtv_filter_analysis(G,g,ft,x,param);

errors = errors + gsp_assert_test(c_exact,c_cheby,1e-4,'JTV - cheby analysis 2');

end

function errors = test_jtv_filter_synthesis()
errors = 0;

N = 500;
T = 300;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

alpha=[0.1:0.1:0.7];
[g,ft] = gsp_jtv_design_wave(G,alpha);

param.lag=1;
x = rand(N,T,numel(alpha));
param.method = 'exact';
s = gsp_jtv_filter_synthesis(G,g,ft,x,param);


param.filtertype = ft;
s2 = gsp_jtv_filter_synthesis_simple(G,g,x,param);

errors = errors + gsp_assert_test(s,s2,1e-10,'JTV  - exact synthesis 1');


param.method = 'cheby';
param.order=200;
s2 = gsp_jtv_filter_synthesis(G,g,ft,x,param);

errors = errors + gsp_assert_test(s,s2,1e-10,'JTV - cheby synthesis 1');

param.method = 'exact';
g = @(x,y) exp(-x.*(y.^2)/0.1);
ft = 'js';
x = rand(N,T,1,3);
s = gsp_jtv_filter_synthesis(G,g,ft,x,param);

param.filtertype = ft;
s2 = gsp_jtv_filter_synthesis_simple(G,g,x,param);
errors = errors + gsp_assert_test(s,s2,1e-10,'JTV  - exact synthesis 2');


param.method = 'cheby';
param.order=20;
s2 = gsp_jtv_filter_synthesis(G,g,ft,x,param);
time=toc;

errors = errors + gsp_assert_test(s,s2,1e-4,'JTV - cheby analysis 2');


end

function errors = test_jtv_matrix_analysis_op()
errors = 0;

N = 20;
T = 20;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

alpha = 1;
[g,ft] = gsp_jtv_design_wave(G,alpha);

F = gsp_jtv_compute_frame(G,g,ft);

x = rand(N,T);

c1 = gsp_jtv_frame_analysis(F,x);

param.method = 'exact';
c2 = gsp_jtv_filter_analysis(G,g,ft,x,param);

errors = errors + gsp_assert_test(c1,c2,eps(1000),'JTV - analysis operator matrix');

end

function errors = test_jtv_matrix_synthesis_op()
errors = 0;

N = 20;
T = 20;
G = gsp_sensor(N);
param.extension=1;
G = gsp_jtv_graph(G,T,[],param);
G = gsp_compute_fourier_basis(G);

[g,ft] = gsp_jtv_design_wave(G);

F = gsp_jtv_compute_frame(G,g,ft);

x = randn(N,T);

param.method = 'exact';

coeff = gsp_jtv_filter_analysis(G,g,ft,x,param);

x1 = gsp_jtv_frame_synthesis(F,coeff);

x2 = gsp_jtv_filter_synthesis(G,g,ft,coeff,param);

errors = errors + gsp_assert_test(x1,x2,eps(1000),'JTV - synthesis operator matrix');

end


function errors = test_jtv_filter_dual()
errors = 0;


N = 50;
T = 50;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

alpha=linspace(0,2,5);
[g,ft] = gsp_jtv_design_damped_wave(G,alpha,[0.01 0.1]);

try
    gd = gsp_jtv_design_can_dual(g,ft);
    gsp_plot_jtv_filter(G,gd,ft);
    close
    fprintf('Test JTV - Canonical 1: OK \n');
catch
    errors = errors +1;
    warning('Test JTV - Error canonical dual')
end


N = 50;
T = 50;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

Nf = 4;
g = gsp_jtv_design_meyer(G,Nf);

try
    gd = gsp_jtv_design_can_dual(g,ft);
    gsp_plot_jtv_filter(G,gd,ft);
    close
    fprintf('Test JTV - Canonical 2: OK \n');
catch
    errors = errors +1;
    warning('Test JTV - Error canonical dual')
end




end



function errors = test_jtv_filter_inverse()
%%
errors = 0;

N = 50;
T = 50;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);

alpha=linspace(0,2,5);
[g,ft] = gsp_jtv_design_damped_wave(G,alpha);
f = randn(N,T);

f1 = gsp_jtv_filter_analysis(G,g,ft,f);
f2 = gsp_jtv_filter_inverse(G,g,ft,f1);

errors = errors + gsp_assert_test(f,f2,eps(1000),'JTV - filter inverse 1');

[h,ft] = gsp_jtv_filter_array(G,g,ft);
f1 = gsp_jtv_filter_analysis(G,h,ft,f);
f2 = gsp_jtv_filter_inverse(G,h,ft,f1);
errors = errors + gsp_assert_test(f,f2,eps(1000),'JTV - filter inverse 2');

Nf = 4;
[g,ft] = gsp_jtv_design_meyer(G,Nf);
f1 = gsp_jtv_filter_analysis(G,g,ft,f);
f2 = gsp_jtv_filter_inverse(G,g,ft,f1);

errors = errors + gsp_assert_test(f,f2,eps(1000),'JTV - filter inverse 3');

[h,ft] = gsp_jtv_filter_array(G,g,ft);
f1 = gsp_jtv_filter_analysis(G,h,ft,f);
f2 = gsp_jtv_filter_inverse(G,h,ft,f1);
errors = errors + gsp_assert_test(f,f2,eps(1000),'JTV - filter inverse 4');

end

function errors = test_theory()


errors = 0;

N = 100;
T = 50;
G = gsp_sensor(N);
G = gsp_jtv_graph(G,T);
G = gsp_compute_fourier_basis(G);
g = @(x,w) abs(exp(-w.*x)).^2;
ft = 'js';
garray = gsp_jtv_filter_array(G,g,ft);
x1p = gsp_itft(G, gsp_vec2mat(gsp_filter_analysis(G,garray,gsp_delta(G,1) ),T))/sqrt(T);
x2p = gsp_jtv_filter_analysis(G,g,ft,gsp_jtv_delta(G,1,1));
% norm(x1p(:)/norm(x1p(:))-x2p(:)/norm(x2p(:)))

errors = errors +gsp_assert_test(x1p,x2p,1e-10,'TEST THEORY')

end
