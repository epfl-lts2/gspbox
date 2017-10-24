function errors = test_gsp_proj_filerbank()

errors = 0;

N = 100;
Nf = 10;

G = gsp_sensor(N);

G = gsp_compute_fourier_basis(G);

g = gsp_design_itersine(G,Nf);

y = rand(N,1);
x = rand(N*Nf,1);

s1 = gsp_proj_filterbank(x,0,G,g,y);
y2 = gsp_filter_synthesis(G,g,s1);
if norm(y-y2)/norm(y)<1e-6
    fprintf('  Test gsp_proj_filterbank OK\n')
else
    fprintf('  Test gsp_proj_filterbank Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end

Nrep = 10;
x = rand(N*Nf,Nrep);

s1 = gsp_proj_filterbank(x,0,G,g,y);

y2 = gsp_filter_synthesis(G,g,s1);

if norm(repmat(y,1,Nrep)-y2,'fro')/norm(y)<1e-6
    fprintf('  Test gsp_proj_filterbank multi OK\n')
else
    fprintf('  Test gsp_proj_filterbank multi Pas OK!!!!!!!!!!!!!!!!\n')
    norm(repmat(y,1,Nrep)-y2,'fro')/norm(y)<1e-6
    errors= errors +1;
end

x = rand(N*Nf,1);
s1 = gsp_proj_filterbank(x,0,G,g,y);
F = gsp_filterbank_matrix(G,g);
param.method = 'exact';
param.A = F;
param.y = y;
param.pinvA = pinv(F);
s2 = proj_linear_eq(x,0,param);

if norm(s1-s2)/norm(s2)<1e-6
    fprintf('  Test gsp_proj_filterbank exact OK\n')
else
    fprintf('  Test gsp_proj_filterbank exact Pas OK!!!!!!!!!!!!!!!!\n')
    norm(s1-s2)/norm(s2)
    errors= errors +1;
end


N = 100;

G = gsp_sensor(N);

G = gsp_compute_fourier_basis(G);

g = gsp_design_expwin(G,0.1);
a = randn(N,1);
y = gsp_filter_analysis(G,g, a);
x = rand(N,1);

s1 = gsp_proj_filterbank(x,0,G,g,y);
y2 = gsp_filter_synthesis(G,g,s1);
if norm(y-y2)/norm(y)<1e-6
    fprintf('  Test gsp_proj_filterbank Fourier 1 OK\n')
else
    fprintf('  Test gsp_proj_filterbank Fourier 1 Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end


N = 1000;

G = gsp_sensor(N);

G = gsp_estimate_lmax(G);

param = struct;
param.method = 'cheby';
param.order = 100;
g = gsp_design_expwin(G,0.1);
g = gsp_design_simple_tf(G,5);
a = randn(N*numel(g),1);
y = gsp_filter_synthesis(G,g, a,param);
x = rand(N*numel(g),1);

s1 = gsp_proj_filterbank(x,0,G,g,y,param);
y2 = gsp_filter_synthesis(G,g,s1,param);
if norm(y-y2)/norm(y)<5e-3
    fprintf('  Test gsp_proj_filterbank cheby 1 OK\n')
else
    fprintf('  Test gsp_proj_filterbank cheby 1 Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end


N = 1000;

G = gsp_sensor(N);

G = gsp_estimate_lmax(G);

param = struct;
param.method = 'lanczos';
param.order = 100;
g = gsp_design_expwin(G,0.1);
g = gsp_design_simple_tf(G,5);
a = randn(N*numel(g),1);
y = gsp_filter_synthesis(G,g, a,param);
x = rand(N*numel(g),1);

s1 = gsp_proj_filterbank(x,0,G,g,y,param);
y2 = gsp_filter_synthesis(G,g,s1,param);
if norm(y-y2)/norm(y)<3e-3
    fprintf('  Test gsp_proj_filterbank lanczos 1 OK\n')
else
    fprintf('  Test gsp_proj_filterbank lanczos 1 Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end



N = 1000;

G = gsp_sensor(N);

G = gsp_estimate_lmax(G);
G = gsp_compute_fourier_basis(G);

param = struct;
param.proj_method = 'proj_b2';
g = gsp_design_expwin(G,0.1);
a = randn(N*numel(g),1);
y = gsp_filter_synthesis(G,g, a);
x = rand(N*numel(g),1);

s1 = gsp_proj_filterbank(x,0,G,g,y,param);
y2 = gsp_filter_synthesis(G,g,s1);
if norm(y-y2)/norm(y)<3e-2
    fprintf('  Test gsp_proj_filterbank proj_b2 1 OK\n')
else
    fprintf('  Test gsp_proj_filterbank proj_b2 1 Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end


N = 1000;

G = gsp_sensor(N);

G = gsp_estimate_lmax(G);

param = struct;
param.proj_method = 'proj_b2';
g = gsp_design_expwin(G,0.1);
a = randn(N*numel(g),1);
y = gsp_filter_synthesis(G,g, a);
x = rand(N*numel(g),1);

s1 = gsp_proj_filterbank(x,0,G,g,y,param);
y2 = gsp_filter_synthesis(G,g,s1);
if norm(y-y2)/norm(y)<3e-2
    fprintf('  Test gsp_proj_filterbank proj_b2 2 OK\n')
else
    fprintf('  Test gsp_proj_filterbank proj_b2 2 Pas OK!!!!!!!!!!!!!!!!\n')
    norm(y-y2)/norm(y)
    errors= errors +1;
end

end