function [ errors ] = test_gsp_filter_manifold(  )

errors = 0;
errors = errors + test_combinatorial();
errors = errors + test_normalized();



end

function errors = test_combinatorial()
%TEST_GSP_FILTER_MANIFOLD 
errors = 0;
gsp_reset_seed(0);
%%
Nx = 20+2;
Ny = Nx;
N = 100;


[X,Y] = meshgrid(linspace(0,1,Nx),linspace(0,1,Ny));
 
X = X(:);
Y = Y(:);

p = [X';Y'];
%%
%G = gsp_2dgrid(sqrt(N));

G = gsp_sensor(N);

t = 0.5/G.N^(1.9/2);


[ G ] = gsp_create_cont_expW( G.coords,G.coords,t );
G = gsp_compute_fourier_basis(G);


G = gsp_estimate_lmax(G);
G.lmax = G.Gm.lmax;
 Nf = 7;
g = gsp_design_mexican_hat(G, Nf); 


sig = rand(G.N,10);
param.method = 'cheby';
r = gsp_filter_analysis(G, g, sig,param);
paramt.method = 'cheby';
rp = gsp_filter_analysis_manifold(G, g,zeros(size(sig)), sig,paramt);


if norm(r-rp,'fro') > 1e-12
    errors = errors + 1;  
    norm(r-rp,'fro')
    fprintf('GSP MANIFOLD ANALYSIS 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD ANALYSIS 1 OK\n')
end


%%


[ G ] = gsp_create_cont_expW( G.coords,G.coords,t );
G = gsp_estimate_lmax(G);
G.lmax = G.Gm.lmax;

g = gsp_design_mexican_hat(G, Nf); 

r3j = rand(G.N,Nf);



r3j = gsp_mat2vec(r3j);

smg = gsp_filter_synthesis(G, g, r3j,param);
%gsp_plot_signal(G,smg)


[smg3, smg2] = gsp_filter_synthesis_manifold(G, g, r3j,r3j);


if norm(smg-smg2)/norm(smg) > 1e-12
    errors = errors + 1;  
     norm(smg-smg2)/norm(smg)
    fprintf('GSP MANIFOLD Syntesis 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis 1 OK\n')
end


if norm(smg-smg3)/norm(smg) > 1e-12
    errors = errors + 1;  
     norm(smg-smg3)/norm(smg)
    fprintf('GSP MANIFOLD Syntesis 2 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis 2 OK\n')
end
%%


N = 100;
Nx = 5;

paramg.distribute = 1;
G = gsp_sensor(N,paramg);
vx = linspace(min(G.coords(:,1)),max(G.coords(:,1)),Nx);
vy = linspace(min(G.coords(:,2)),max(G.coords(:,2)),Nx);
[X, Y] = meshgrid(vx, vy );
X = X(:);
Y = Y(:);
p = [X';Y'];
%p = [0.5; 0.5];
t = 0.005;
[ G ] = gsp_create_cont_expW(G.coords, p' ,t);

G = gsp_compute_fourier_basis(G);


g = gsp_design_heat(G,50);

% To compute all atoms
paramt.method = 'cheby';
M11 = gsp_filter_analysis_manifold(G, g,zeros(size(p,2),G.N),eye(G.N),paramt);
[~,M21] = gsp_filter_synthesis_manifold(G, g,zeros(size(p,2),G.N),eye(G.N),paramt);


if norm(M11-M21,'fro') > 1e-12
    errors = errors + 1;  
     norm(M11-M21,'fro')
    fprintf('GSP MANIFOLD sym 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD sym 1 OK\n')
end

paramt.method = 'cheby';
xg = rand(G.N,3);
xout = rand(size(p,2),3);
[c11,c12] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);
paramt.method = 'exact';
[c21,c22] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);


c1 = [c11;c12];

c2 = [c21;c22];

G = gsp_build_oose_fourier_basis(G);


if norm(G.Gm.L*G.Gm.U-(diag(G.Gm.e)*G.Gm.U')','fro')/norm(G.Gm.L,'fro') > 1e-10
    errors = errors + 1;  
    norm(G.Gm.L*G.Gm.U-(diag(G.Gm.e)*G.Gm.U')','fro')/norm(G.Gm.L,'fro')
    fprintf('GSP MANIFOLD Fourier basis PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis OK\n')
end

if norm(G.Gm.U*diag(G.Gm.e)*G.Gm.Um1 - G.Gm.L,'fro')/norm(G.Gm.L,'fro') > 1e-7
    errors = errors + 1;  
    norm(G.Gm.U*diag(G.Gm.e)*G.Gm.Um1 - G.Gm.L,'fro')/norm(G.Gm.L,'fro')
    fprintf('GSP MANIFOLD Fourier basis 2 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis 2 OK\n')
end

if norm(c1-c2)/norm(c1) > 1e-7
    errors = errors + 1;  
     norm(c1-c2,'fro')/norm(c1)
    fprintf('GSP MANIFOLD exact 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD exact 1 OK\n')
end

% paramt.method = 'lanczos';
% [c31,c32] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);
% c3 = [c31;c32];
% 
% if norm(c1-c3)/norm(c1) > 1e-7
%     errors = errors + 1;  
%      norm(c1-c3,'fro')/norm(c1)
%     fprintf('GSP MANIFOLD lanczos 1 PAS OK!!!!!!!!!!!!!!!!!\n')
% else
%     fprintf('GSP MANIFOLD lanczos 1 OK\n')
% end
% 


paramt.method = 'cheby';
paramt.order =100;
xg = rand(G.N*length(g),5);
xout = zeros(size(p,2)*length(g),5);
[c11,c12] = gsp_filter_synthesis_manifold(G, g,xout,xg,paramt);
paramt.method = 'exact';
[c21,c22] = gsp_filter_synthesis_manifold(G, g,xout,xg,paramt);


c1 = [c11;c12];

c2 = [c21;c22];


if norm(c1-c2)/norm(c1) > 1e-8
    errors = errors + 1;  
     norm(c1-c2)/norm(c1)
    fprintf('GSP MANIFOLD Syntesis exact PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis exact OK\n')
end




paramt.method = 'exact';
paramt.order =30;
xg = rand(G.N*length(g),5);
xout = rand(size(p,2)*length(g),5);
%xout = zeros(size(p,2)*length(g),5);


g = gsp_design_mexican_hat(G,7);

c1 = gsp_filter_analysis_manifold(G, g,xout, xg,paramt);

Fr = gsp_filterbank_matrix_manifold(G,g,paramt);
xtot = [xg;xout];




c2 = Fr' * xtot;
%c2 = gsp_filter_analysis(G,g,xg,paramt);

if norm(c1-c2)/norm(c1) > 1e-8
    errors = errors + 1;  
     norm(c1-c2)/norm(c1)
    fprintf('GSP MANIFOLD filterbank matrix PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD filterbank matrix OK\n')
end


end




function errors = test_normalized()
%TEST_GSP_FILTER_MANIFOLD 
errors = 0;
gsp_reset_seed(0);
%%
Nx = 20+2;
Ny = Nx;
N = 100;


[X,Y] = meshgrid(linspace(0,1,Nx),linspace(0,1,Ny));
 
X = X(:);
Y = Y(:);

%p = [X';Y'];
%%
G = gsp_2dgrid(sqrt(N));

G = gsp_sensor(N);

t = 0.5/G.N^(1.9/2);


[ G ] = gsp_create_cont_expW( G.coords,G.coords,t );
G = gsp_create_laplacian(G,'normalized');
G = gsp_compute_fourier_basis(G);


 Nf = 7;
g = gsp_design_mexican_hat(G, Nf); 


sig = rand(G.N,10);
param.method = 'cheby';
r = gsp_filter_analysis(G, g, sig,param);
paramt.method = 'cheby';
rp = gsp_filter_analysis_manifold(G, g,zeros(size(sig)), sig,paramt);


if norm(r-rp,'fro') > 1e-12
    errors = errors + 1;  
    norm(r-rp,'fro')
    fprintf('GSP MANIFOLD ANALYSIS NORMALIZED 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD ANALYSIS NORMALIZED 1 OK\n')
end


%%


[ G ] = gsp_create_cont_expW( G.coords,G.coords,t );
G = gsp_create_laplacian(G,'normalized');
G = gsp_estimate_lmax(G);

g = gsp_design_mexican_hat(G, Nf); 

r3j = rand(G.N,Nf);



r3j = gsp_mat2vec(r3j);

smg = gsp_filter_synthesis(G, g, r3j,param);
%gsp_plot_signal(G,smg)


[smg3, smg2] = gsp_filter_synthesis_manifold(G, g, r3j,r3j);


if norm(smg-smg2)/norm(smg) > 1e-12
    errors = errors + 1;  
     norm(smg-smg2)/norm(smg)
    fprintf('GSP MANIFOLD Syntesis 1 NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis 1 NORMALIZED OK\n')
end


if norm(smg-smg3)/norm(smg) > 1e-12
    errors = errors + 1;  
     norm(smg-smg3)/norm(smg)
    fprintf('GSP MANIFOLD Syntesis 2 NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis 2 NORMALIZED OK\n')
end
%%


N = 100;
Nx = 5;

paramg.distribute = 1;
G = gsp_sensor(N,paramg);
vx = linspace(min(G.coords(:,1)),max(G.coords(:,1)),Nx);
vy = linspace(min(G.coords(:,2)),max(G.coords(:,2)),Nx);
[X, Y] = meshgrid(vx, vy );
X = X(:);
Y = Y(:);
p = [X';Y'];
%p = [0.5; 0.5];
t = 0.005;
[ G ] = gsp_create_cont_expW( G.coords,p',t);
[ G ] = gsp_create_laplacian(G,'normalized');

G = gsp_compute_fourier_basis(G);


g = gsp_design_heat(G,50);

% To compute all atoms
paramt.method = 'cheby';
M11 = gsp_filter_analysis_manifold(G, g,zeros(size(p,2),G.N),eye(G.N),paramt);
[~,M21] = gsp_filter_synthesis_manifold(G, g,zeros(size(p,2),G.N),eye(G.N),paramt);


if norm(M11-M21,'fro') > 1e-12
    errors = errors + 1;  
     norm(M11-M21,'fro')
    fprintf('GSP MANIFOLD sym 1 NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD sym 1 NORMALIZED OK\n')
end

paramt.method = 'cheby';
xg = rand(G.N,3);
xout = rand(size(p,2),3);
[c11,c12] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);
paramt.method = 'exact';
[c21,c22] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);


c1 = [c11;c12];

c2 = [c21;c22];
G2 = G;
G = gsp_build_oose_fourier_basis(G);

G2 = gsp_build_oose_fourier_basis(G2,1);



if norm(G.Gm.L*G.Gm.U-(diag(G.Gm.e)*G.Gm.U')','fro')/norm(G.Gm.L,'fro') > 1e-10
    errors = errors + 1;  
    norm(G.Gm.L*G.Gm.U-(diag(G.Gm.e)*G.Gm.U')','fro')/norm(G.Gm.L,'fro')
    fprintf('GSP MANIFOLD Fourier basis PAS NORMALIZED OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis OK\n')
end

if norm(G.Gm.U*diag(G.Gm.e)*G.Gm.Um1 - G.Gm.L,'fro')/norm(G.Gm.L,'fro') > 1e-7
    errors = errors + 1;  
    norm(G.Gm.U*diag(G.Gm.e)*G.Gm.Um1 - G.Gm.L,'fro')/norm(G.Gm.L,'fro')
    fprintf('GSP MANIFOLD Fourier basis 2 NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis 2 NORMALIZED OK\n')
end

if norm(G2.Gm.U*diag(G2.Gm.e)*G2.Gm.Um1 - G2.Gm.L,'fro')/norm(G2.Gm.L,'fro') > 1e-7
    errors = errors + 1;  
    norm(G2.Gm.U*diag(G2.Gm.e)*G2.Gm.Um1 - G2.Gm.L,'fro')/norm(G2.Gm.L,'fro')
    fprintf('GSP MANIFOLD Fourier basis sparse 1 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis sparse 1 OK\n')
end

if norm(G2.Gm.U - G.Gm.U,'fro')> 1e-7
    errors = errors + 1;  
    fprintf('GSP MANIFOLD Fourier basis sparse 2 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis sparse 2 OK\n')
end

if norm(G2.Gm.Um1 - G.Gm.Um1,'fro')> 1e-7
    errors = errors + 1;  
    fprintf('GSP MANIFOLD Fourier basis sparse 3 PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Fourier basis sparse 3 OK\n')
end

if norm(c1-c2)/norm(c1) > 1e-7
    errors = errors + 1;  
     norm(c1-c2,'fro')/norm(c1)
    fprintf('GSP MANIFOLD exact 1 NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD exact 1 NORMALIZED OK\n')
end

% paramt.method = 'lanczos';
% [c31,c32] = gsp_filter_analysis_manifold(G, g,xg,xout,paramt);
% c3 = [c31;c32];
% 
% if norm(c1-c3)/norm(c1) > 1e-7
%     errors = errors + 1;  
%      norm(c1-c3,'fro')/norm(c1)
%     fprintf('GSP MANIFOLD lanczos 1 PAS OK!!!!!!!!!!!!!!!!!\n')
% else
%     fprintf('GSP MANIFOLD lanczos 1 OK\n')
% end
% 


paramt.method = 'cheby';
paramt.order =100;
xg = rand(G.N*length(g),5);
xout = zeros(size(p,2)*length(g),5);
[c11,c12] = gsp_filter_synthesis_manifold(G, g,xout, xg,paramt);
paramt.method = 'exact';
[c21,c22] = gsp_filter_synthesis_manifold(G, g,xout,xg,paramt);


c1 = [c11;c12];

c2 = [c21;c22];


if norm(c1-c2)/norm(c1) > 1e-8
    errors = errors + 1;  
     norm(c1-c2)/norm(c1)
    fprintf('GSP MANIFOLD Syntesis exact NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD Syntesis exact NORMALIZED OK\n')
end




paramt.method = 'exact';
paramt.order =30;
xg = rand(G.N*length(g),5);
xout = rand(size(p,2)*length(g),5);
%xout = zeros(size(p,2)*length(g),5);


g = gsp_design_mexican_hat(G,7);

c1 = gsp_filter_analysis_manifold(G, g,xout, xg,paramt);

Fr = gsp_filterbank_matrix_manifold(G,g,paramt);
xtot = [xg;xout];




c2 = Fr' * xtot;
%c2 = gsp_filter_analysis(G,g,xg,paramt);

if norm(c1-c2)/norm(c1) > 1e-8
    errors = errors + 1;  
     norm(c1-c2)/norm(c1)
    fprintf('GSP MANIFOLD filterbank matrix NORMALIZED PAS OK!!!!!!!!!!!!!!!!!\n')
else
    fprintf('GSP MANIFOLD filterbank matrix NORMALIZED OK\n')
end


end
