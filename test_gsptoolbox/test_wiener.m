
function errors = test_wiener()

errors = 0;
gsp_reset_seed
N = 20;
M = 3;
sigma = 0.3;

G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
g = @(x) sin(x);
psd = @(x) g(x).^2;

Mask = rand(N,1)>0.5;

x = gsp_filter_analysis(G,g,randn(N,M));

Mop =@(x) bsxfun(@times,Mask,x);

y1 = Mop(x);

y2 = Mop(x+sigma*randn(N,M));
clear paramopt
paramopt.maxit = 1000;
paramopt.tol = 1e-12;
sol11 = gsp_wiener_inpainting(G,y1,Mask,psd,0,paramopt);
sol12 = gsp_wiener_inpainting_exact(G,y1,Mask,psd,0);

if norm(sol11-sol12,'fro')/norm(sol11,'fro') > 1e-10
    errors = errors + 1;
    norm(sol11-sol12,'fro')/norm(sol11,'fro') 
    warning('Test Wiener opt 1 - ERROR');
else
    disp('Test Wiener opt 1 - OK');
end


sol21 = gsp_wiener_inpainting(G,y2,Mask,psd,sigma.^2,paramopt);
sol22 = gsp_wiener_inpainting_exact(G,y2,Mask,psd,sigma.^2);

if norm(sol21-sol22,'fro')/norm(sol22,'fro')> 1e-10
    errors = errors + 1;
    norm(sol21-sol22,'fro')/norm(sol22,'fro')
    warning('Test Wiener opt 2 - ERROR');
else
    disp('Test Wiener opt 2 - OK');
end




%%

g = gsp_design_expwin(G,0.2);
psd = @(x) g(x).^2;

Mask = rand(N,1)>0.5;
x = gsp_filter_analysis(G,g,randn(N,M));

Mop =@(x) bsxfun(@times,Mask,x);

y1 = Mop(x);

y2 = Mop(x+sigma*randn(N,M));

paramopt.tol = 1e-12;

%%
paramopt.maxit = 1000;
paramopt.gamma = 0.01;
sol11 = gsp_wiener_inpainting(G,y1,Mask,psd,0,paramopt);
sol12 = gsp_wiener_inpainting_exact(G,y1,Mask,psd,0);

if norm(sol11-sol12,'fro')/norm(sol11,'fro') > 1e-3
    errors = errors + 1;
    norm(sol11-sol12,'fro')/norm(sol11,'fro') 
    warning('Test Wiener opt 3 - ERROR');
else
    disp('Test Wiener opt 3 - OK');
end

if norm(x-sol12,'fro')/norm(x,'fro') > 1e-10
    errors = errors + 1;
    norm(x-sol12,'fro')/norm(x,'fro') 
    warning('Test Wiener opt 4 - ERROR');
else
    disp('Test Wiener opt 4 - OK');
end




%%
% figure(1)
% subplot(121)
% gsp_plot_signal_spectral(G,gsp_gft(G,x(:,1)))
% subplot(122)
% gsp_plot_signal_spectral(G,gsp_gft(G,sol12(:,1)))


%%
paramopt.gamma = 0.1;

paramopt.maxit = 1000;

sol21 = gsp_wiener_inpainting(G,y2,Mask,psd,sigma.^2,paramopt);
sol22 = gsp_wiener_inpainting_exact(G,y2,Mask,psd,sigma.^2);


if norm(sol21-sol22,'fro')/norm(sol22,'fro')> 1e-10
    errors = errors + 1;
    norm(sol21-sol22,'fro')/norm(sol22,'fro')
    warning('Test Wiener opt 5 - ERROR');
else
    disp('Test Wiener opt 5 - OK');
end

end