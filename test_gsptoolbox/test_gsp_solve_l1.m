function [ errors ] = test_gsp_solve_l1(  )

% Testin of gsp_solve_l1
gsp_reset_seed(0)
errors = 0;

N = 100;
Nf = 5;
G = gsp_sensor(N);
%G = gsp_2dgrid(sqrt(N));
G = gsp_compute_fourier_basis(G);

g = gsp_design_itersine(G,Nf);

s = gsp_filter_analysis(G,g{1},randn(N,1));
[alpha1, info1]  = gsp_solve_l1(G, g, s, 0 );

sp = gsp_filter_synthesis(G,g,alpha1);


if norm(sp-s)/norm(s)<1e-8
    fprintf('Test SOLVE_L1 OK\n')
else
    errors = errors +1;
    fprintf('Test SOLVE_L1 PAS OK\n')
    norm(sp-s)/norm(s)
end


[alpha2, info2]  = gsp_solve_l1(G, g, s, 0.01 );
sp = gsp_filter_synthesis(G,g,alpha2);
if norm(alpha2,1)<norm(alpha1,1)
    fprintf('Test SOLVE_L1 2 OK\n')
else
    errors = errors +1;
    fprintf('Test SOLVE_L1 2 PAS OK\n')
    norm(alpha2,1)
    norm(alpha1,1)
end

sp2 = gsp_filter_synthesis(G,g,alpha2);
if norm(sp2-s)/norm(s)<0.05
    fprintf('Test SOLVE_L1 3 OK\n')
else
    errors = errors +1;
    fprintf('Test SOLVE_L1 3 PAS OK\n')
    norm(sp2-s)/norm(s)
end


% old code... probably stupid
% % Testin of gsp_solve_l1
% gsp_reset_seed(0)
% errors = 0;
% 
% N = 100;
% Nf = 5;
% param = struct;
% %param.distribute = 1;
% G = gsp_sensor(N,param);
% %G = gsp_2dgrid(sqrt(N));
% G = gsp_compute_fourier_basis(G);
% G = gsp_spectrum_cdf_approx(G);
% paramfilter.filter = gsp_design_meyer(G,Nf);
% g = gsp_design_warped_translates(G, Nf, paramfilter);   
% %g = gsp_design_itersine(G,Nf);
% 
% % figure
% % gsp_plot_filter(G,g)
% 
% alphain = zeros(N,Nf);
% 
% for ii = 1:4
%     alphain(ii,ii) = 1;
% end
% 
% alphain = gsp_mat2vec(alphain);
% signal = gsp_filter_synthesis(G,g,alphain);
% 
% % figure
% % gsp_plot_signal(G, signal);
% param.verbose = 0;
% param.maxit = 1000;
% % fig=figure(100);
% % param.do_sol=@(x) plot_objective(x,fig);
% 
% param.tol = 1e-6;
% param.guess = alphain;
% lambda = 0.05;
% alphaout = gsp_solve_l1(G,g,signal,lambda,param);
% 
% % param.guess = alphaout;
% % lambda2 = 0.1;
% % param.tol = 1e-10;
% % param.maxit = 200;
% % alphaout = gsp_solve_l0(G,g,signal,lambda2,param);
% % norm(alphaout-alphain)
% 
% if norm(alphaout-alphain)<2
%     fprintf('Test SOLVE_L1 OK\n')
% else
%     errors = errors +1;
%     fprintf('Test SOLVE_L1 PAS OK\n')
%     norm(alphaout-alphain)
% end
    



end

