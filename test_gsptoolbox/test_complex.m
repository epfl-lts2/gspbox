function err = test_complex()

err = 0;

gsp_reset_seed(0)
N = 64;
M = 2;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
% Let us select a complex kernel
k = @(x) 1./(1i+x);
% Let us create the dual kernel
%kd = @(x) conj((k(x).^(-1))); 
kd = gsp_design_can_dual(k);
kabs = @(x) abs(k(x));
kdabs = gsp_design_can_dual(kabs);



%% 1) try on the normal graph

% Let us select a complex signal
x = randn(N,M) + 1i* randn(N,M);
param.method = 'exact';
c = gsp_filter_analysis(G,k,x,param);


xt = gsp_filter_synthesis(G,kd,c,param);

err = err + double(norm(x-xt,'fro')>eps(1e5));
% norm(x-xt,'fro')

G2 = gsp_ring(N);
G2 = gsp_compute_fourier_basis(G2);

%% 2) test analysis on a ring

x = randn(N,M) + 1i* randn(N,M);
param.method = 'exact';
c3 = gsp_filter_analysis(G2,k,x,param);
v = gsp_classic2graph_eig_order( N );
fit = zeros(N,1);
fit(v) = conj(gsp_filter_evaluate(k,G2.e));
c3p = ifft(fft(x).*repmat(fit,1,M));
err = err + double(norm(c3-c3p,'fro')>eps(1e5));

% norm(c3-c3p,'fro')
err = err + double(norm(c3-c3p,'fro')>eps(1e5));

%% 3) test synthesis on a ring
param.method = 'exact';
Amat = gsp_filter_analysis(G2,k,eye(N),param);
param.method = 'exact';
s1 = gsp_filter_synthesis(G2,k,c3,param);
s2 = Amat'*c3;
% norm(s1-s2,'fro')
err = err + double(norm(s1-s2,'fro')>eps(1e5));



%% 4) Try on a ring with real filter
x = randn(N,M) + 1i* randn(N,M);
param.method = 'exact';
c1 = gsp_filter_analysis(G2,kabs,x,param);
xt1 = gsp_filter_synthesis(G2,kdabs,c1,param);

err = err + double(norm(x-xt1,'fro')>eps(1e5));
% norm(x-xt1,'fro')

%% 5) Try on a ring with real signal
x = randn(N,M);
param.method = 'exact';
c2 = gsp_filter_analysis(G2,k,x,param);
xt2 = gsp_filter_synthesis(G2,kd,c2,param);

err = err + double(norm(x-xt2,'fro')>eps(1e5));

% norm(x-xt2,'fro')

%% 5) Try on a ring with everything complex
x = randn(N,M) + 1i* randn(N,M);
param.method = 'exact';
c2 = gsp_filter_analysis(G2,k,x,param);
xt2 = gsp_filter_synthesis(G2,kd,c2,param);

err = err + double(norm(x-xt2,'fro')>eps(1e5));

% norm(x-xt2,'fro')

%% 7) Test pseudo inverse frame

param.method = 'exact';
Amat = gsp_filter_analysis(G2,k,eye(N),param);
param.method = 'exact';
Amatp = gsp_filter_synthesis(G2,kd,eye(N*numel(k)),param);
% norm(pinv(Amat)-Amatp,'fro')
err = err + double(norm(pinv(Amat)-Amatp,'fro')>eps(1e5));

% %% 8) Test what happen with chebyshef and lanczos
% %   real kernel, complexe signal
% % => result lanczos fail!
% 
% x = randn(N,M) + 1i* randn(N,M);
% 
% 
% param.method = 'exact';
% c1 = gsp_filter_analysis(G2,kabs,x,param);
% param.method = 'cheby';
% param.order = 1000;
% c2 = gsp_filter_analysis(G2,kabs,x,param);
% param.method = 'lanczos';
% param.order = 1000;
% c3 = gsp_filter_analysis(G2,kabs,x,param);
% norm(c1-c2,'fro')
% norm(c1-c3,'fro')
% norm(c2-c3,'fro')
% err = err + double(norm(c-c2,'fro')>eps(1e5));
% err = err + double(norm(c1-c3,'fro')>eps(1e5));

% %% 8) Test what happen with chebyshef
% % Real signal, complex kernel
% % Result => both fail!
% param.method = 'exact';
% Amat = gsp_filter_analysis(G2,k,eye(N),param);
% param.method = 'cheby';
% param.order = 1000;
% Amat2 = gsp_filter_analysis(G2,k,eye(N),param);
% param.method = 'lanczos';
% param.order = 1000;
% Amat3 = gsp_filter_analysis(G2,k,eye(N),param);
% norm(Amat-Amat2,'fro')
% norm(Amat-Amat3,'fro')
% norm(Amat2-Amat3,'fro')
% err = err + double(norm(Amat-Amat2,'fro')>eps(1e5));
% err = err + double(norm(Amat-Amat3,'fro')>eps(1e5));


if err
    warning('error in test_complex')
end

end
