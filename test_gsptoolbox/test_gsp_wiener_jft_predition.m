function errors = test_gsp_wiener_jft_predition
errors = 0;

%% 
%Time parameter
T = 200;
% Graph parameters
N = 100;

param.accelerate = 0;

%Graph
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

G = gsp_jtv_graph(G,T);


%wave filter
alpha=0.01;
g = @(lambda,f) exp(-(lambda-4*G.lmax*f.^2).^2/alpha);

psdfilt = @(lambda,f) abs(g(lambda,f)).^2;

ft = 'js';
psdcoeff = abs(gsp_jtv_filter_evaluate(g,ft,G.e,gsp_jtv_fa(G))).^2;


[psdfiltarray, ftarray] = gsp_jtv_filter_array(G,psdfilt,ft);

%% build the problem

x0 = full(sprandn(N,T,0.01));

X = gsp_jtv_filter_synthesis(G,g,ft,x0);

%%

%M = rand(N,T)>0.5;
M = diag(rand(N,1)>0.5)*ones(N,T);

sigma = 0;
psdnoise = sigma.^2;
y = M.*X + sigma*randn(N,T);
param.maxit = 100;
param.tol = 1e-16;

X1 = gsp_jtv_wiener_inpainting(G,y,M,psdcoeff,psdnoise,param);

X2 = gsp_jtv_wiener_inpainting(G,y,M,psdfilt,psdnoise,param);
errors = errors + gsp_assert_test(X1,X2,1e-10,'WIENER JTV coeff vs filt 1');

%%

%M = rand(N,T)>0.5;
M = diag(rand(N,1)>0.5)*ones(N,T);

sigma = 0.5;
psdnoise = sigma.^2;
y = M.*X + sigma*randn(N,T);

X1 = gsp_jtv_wiener_inpainting(G,y,M,psdcoeff,psdnoise,param);

X2 = gsp_jtv_wiener_inpainting(G,y,M,psdfilt,psdnoise,param);
X3 = gsp_jtv_wiener_inpainting(G,y,M,psdfilt,@(x,T) psdnoise,param);

errors = errors + gsp_assert_test(X1,X2,1e-10,'WIENER JTV coeff vs filt 2');
errors = errors + gsp_assert_test(X1,X3,1e-10,'WIENER JTV coeff vs filt 2 prime');



%%

%M = rand(N,T)>0.5;
M = diag(rand(N,1)>0.5)*ones(N,T);

sigma = 0;
psdnoise = sigma.^2;
y = M.*X + sigma*randn(N,T);
param.maxit = 100;
param.tol = 1e-16;

tmp1 = gsp_jtv_filter_evaluate(psdfiltarray,ftarray,G.e,gsp_jtv_fa(G));
tmp2 = gsp_jtv_filter_evaluate(psdfilt,ft,G.e,gsp_jtv_fa(G));
errors = errors + gsp_assert_test(tmp2,tmp1,1e-10,'WIENER JTV coeff vs filt 3 a');

errors = errors + gsp_assert_test(psdcoeff,tmp2,1e-10,'WIENER JTV coeff vs filt 3 a prime');

fprox1 = @(T) @(lambda,omega) psdfilt(lambda,omega)./(2*T + psdfilt(lambda,omega) + eps);
fprox2= @(T) apply2array(psdfiltarray, @(x) x./(2*T + x + eps));

tmp = rand;
fv1 = fprox1(tmp);
fv2 = fprox2(tmp);



tmp1 = gsp_jtv_filter_evaluate(fv1,ft,G.e,gsp_jtv_fa(G));
tmp2 = gsp_jtv_filter_evaluate(fv2,[ft,'-array'],G.e,gsp_jtv_fa(G));
errors = errors + gsp_assert_test(tmp1,tmp2,1e-10,'WIENER JTV coeff vs filt 3 b');

wl1 = @(lambda,omega) 1./(psdfilt(lambda,omega)+eps);

wl2 = apply2array(psdfiltarray, @(x) 1./(x+eps));

tmp1 = gsp_jtv_filter_evaluate(wl1,ft,G.e,gsp_jtv_fa(G));
tmp2 = gsp_jtv_filter_evaluate(wl2,[ft,'-array'],G.e,gsp_jtv_fa(G));
errors = errors + gsp_assert_test(tmp1,tmp2,1e-10,'WIENER JTV coeff vs filt 3 c');
%%
X1 = gsp_jtv_wiener_inpainting(G,y,M,psdcoeff,psdnoise,param);
X2 = gsp_jtv_wiener_inpainting(G,y,M,psdfiltarray,psdnoise,param);
errors = errors + gsp_assert_test(X1,X2,1e-10,'WIENER JTV coeff vs filt 3 final');

%%

%M = rand(N,T)>0.5;
M = diag(rand(N,1)>0.5)*ones(N,T);

sigma = 0.5;
psdnoise = sigma.^2;
y = M.*X + sigma*randn(N,T);

X1 = gsp_jtv_wiener_inpainting(G,y,M,psdcoeff,psdnoise,param);

X2 = gsp_jtv_wiener_inpainting(G,y,M,psdfiltarray,psdnoise,param);

errors = errors + gsp_assert_test(X1,X2,1e-10,'WIENER JTV coeff vs filt 4');

end




function g = apply2array(garray,fun)

g = cell(size(garray));
for ii = 1:numel(garray)
    g{ii} = @(x) fun(garray{ii}(x));
end

end