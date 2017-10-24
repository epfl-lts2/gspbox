%GSP_DEMO_POSTER A 6 steps demonstration
%    
%   This file contains a 6 steps demonstration of the GSPBOX
%

N = 100; % number of nodes
% 1) Create a graph
G = gsp_sensor(N);
% 2) Compute the Fourier basis
G = gsp_compute_fourier_basis(G);
% 3) Create a smooth signal with noise
x = G.U(:,2);
y = x + 1/sqrt(N)*randn(N,1);
% 4) Select a filter 
g = gsp_design_expwin(G,0.1);
% 5) Remove the noise
s = gsp_filter(G,g,y);
% 6) Display the results
figure(1); gsp_plot_signal(G,x); title('Original signal');
figure(2); gsp_plot_signal(G,y); title('Noisy signal');
figure(3); gsp_plot_signal(G,s); title('Denoised signal');

% Show the filter
figure(4); gsp_plot_filter(G,g)

% Save the plots
figure(1); gsp_plotfig('poster_ori');
figure(2); gsp_plotfig('poster_noisy');
figure(3); gsp_plotfig('poster_denoised');