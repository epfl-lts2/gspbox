function [ errors ] = test_gsp_dn( )

errors = 0;

errors = errors + test_wavelet_dn();



end


function errors = test_wavelet_dn()
errors = 0;

lambda = 1;


Nf = 6;
N = 100;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G); 

W = gsp_design_mexican_hat(G,Nf); % Generate tight frame wavelet

f = sign(G.U(:,2));
sigma = 0.3;
gsp_reset_seed(0);
x = f+sigma*randn(size(f));
param.tol = 1e-12;
param.maxit = 100;
param.verbose = 0;
param.method = 'FISTA';
y1 = gsp_wavelet_dn(G,W,x,lambda,param);
param.method = 'DG';
param.gamma = 1;
y2 = gsp_wavelet_dn(G,W,x,lambda,param);


if norm(y1(:)-y2(:))/norm(y1(:)) > 1e-9
   errors = errors +1;
   warning('DENOISING: Wavlet DG diff FISTA error')
   badness = norm(y1(:)-y2(:))/norm(y1(:))
   figure;
    gsp_plot_signal(G,x);
    figure;
    gsp_plot_signal(G,y1);

else
    fprintf('DENOISING: DG == FISTA ok\n');
end


end
