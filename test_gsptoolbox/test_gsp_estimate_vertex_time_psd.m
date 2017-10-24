function errors = test_gsp_estimate_vertex_time_psd
    
errors = 0;

%Time parameter
T = 200;fs=1;
% Graph parameters
N = 100;
%Graph
G = gsp_sensor(N);
% G1 = gsp_2dgrid(10);
G = gsp_jtv_graph(G,T,fs);
G = gsp_compute_fourier_basis(G);


%wave filter
alpha=1;
beta = 0.7;
[g, ft] = gsp_jtv_design_damped_wave(G, alpha,beta);


Nx = 5;
x2 = (rand(N,T,1,Nx)>0.99).*randn(N,T,1,Nx);
X2 = gsp_jtv_filter_synthesis(G,g,ft,x2);

param.L = 100;
param.use_fast = 0;
t1 = tic;
param.estimator = 'TVA';
psd_TVA1 = gsp_estimate_vertex_time_psd(G,X2,param);
time1 = toc(t1)

t2 = tic;
param.estimator = 'TVA';
param.use_fast = 1;
psd_TVA2 = gsp_estimate_vertex_time_psd(G,X2,param);
time2 = toc(t2)

A1 = gsp_filter_evaluate(psd_TVA1,G.e);
A2 = gsp_filter_evaluate(psd_TVA2,G.e);

gsp_assert_test(A1,A2,1e-10,'ESTIMATE_PSD: TVA')



end