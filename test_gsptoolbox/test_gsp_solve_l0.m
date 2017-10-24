function [ errors ] = test_gsp_solve_l0(  )
% Testin of gsp_solve_l1
errors = 0;

N = 100;
Nf = 5;
param.distribute = 1;
G = gsp_sensor(N,param);
%G = gsp_2dgrid(sqrt(N));
G = gsp_compute_fourier_basis(G);
G = gsp_spectrum_cdf_approx(G);
paramfilter.filter = gsp_design_meyer(G,Nf);
g = gsp_design_warped_translates(G, Nf, paramfilter);   
%g = gsp_design_itersine(G,Nf);

% figure
% gsp_plot_filter(G,g)

alphain = zeros(N,Nf);

for ii = 2:5
    alphain(ii,ii) = 1;
end

alphain = gsp_mat2vec(alphain);
signal = gsp_filter_synthesis(G,g,alphain);

% figure
% gsp_plot_signal(G, signal);
param.verbose = 0;
% fig=figure(100);
% param.do_sol=@(x) plot_objective(x,fig);

lambda = 0.1;

param.tol = 1e-10;
%alphaout1 = gsp_solve_l0(G,g,signal,lambda,param);
param.guess = alphain;
alphaout = gsp_solve_l0(G,g,signal,lambda,param);


if norm(alphaout-alphain)<1e-5
    fprintf('Test SOLVE_L0 OK\n')
else
    errors = errors +1;
    fprintf('Test SOLVE_L0 PAS OK\n')
    norm(alphaout-alphain)
end
    
lambda = 10;
alphaout = gsp_solve_l0(G,g,signal,lambda,param);


if sum(alphaout) == 0
    fprintf('Test SOLVE_L0 2 OK\n')
else
    errors = errors +1;
    fprintf('Test SOLVE_L0 2 PAS OK\n')
    norm(alphaout-alphain)
end

end

