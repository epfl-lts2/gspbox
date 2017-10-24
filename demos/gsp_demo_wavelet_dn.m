%GSP_DEMO_WAVELET_DN Demonstratration of the use of wavelet for denoising
%
%   In this small example, we show how to perform wavelet denoising using
%   the GSPBox and particularly the function |gsp_wavelet_dn| .
%
%   The function |gsp_wavelet_dn| removes the low frequency part of the
%   signal and solves a l1 minimization problem to remove the noise.
%
%   .. figure::
%
%      Result of filtering
%
%      We observe that the wavelet denoising allows rapid change on the
%      graph signal. This is not possible with a simple low pass filtering.
%
%   .. figure::
%
%      Choosen wavelet filterbank
%
%      Here we use a mexican hat frame for the wavelet construction. 
%   

%% Initialization
clear
close all

gsp_reset_seed(0);

%% Parameters

N = 225;           % Number of samples
sigma = 0.25;      % Noise level
Nf = 7;            % Number of filter
Tau = sigma;       % Regularization parameter

% Graph
param_sensor.distribute = 1;
G = gsp_sensor(N,param_sensor);
G.plotting.edge_color = [0.6, 0.6, 0.6];
G = gsp_estimate_lmax(G);
G = gsp_compute_fourier_basis(G);

%% Create a signal

% % Piecewise smooth signal
% f = rand(N,1);
% f = f > (max(f)/(1+1/10));
% g = gsp_design_heat(G,30);
% f = gsp_filter(G,g,f);
% f(f>0.1) = f(f>0.1)+0.1;
% f(f<0.05) = f(f<0.05)-0.1;
% f = f/norm(f)*sqrt(N);

% Cut signal
f = sign(G.U(:,2));

% Add the noise
fn = f + sigma*randn(N,1);


%% Select a filter bank

paramw.normalize = 0;
w = gsp_design_mexican_hat(G,Nf,paramw);

%% Solve the problem
param_dn.maxit = 100;

% When Tau is big, the algorithm just perform a low pass filtering
fdn_lp = gsp_wavelet_dn(G,w,fn,100*Tau,param_dn);
fdn = gsp_wavelet_dn(G,w,fn,Tau,param_dn);

%% Plot the result

figure(1)
subplot(221)
gsp_plot_signal(G,f);
title('Original signal')

subplot(222)
gsp_plot_signal(G,fn);
title('Noisy signal')

subplot(223)
gsp_plot_signal(G,fdn_lp);
title('Denoised signal LP')

subplot(224)
gsp_plot_signal(G,fdn);
title('Denoised signal wavelets')

figure;
gsp_plot_filter(G,w);
title('Wavelet filterbank')
