%GSP_DEMO_WAVELET_DN Demonstratration of the use of wavelet for denoising
%
%   In this small example, we show how to perform wavelet denoising using
%   the GSPBox and particularly the function GSP_WAVELET_DN .
%
%   The function GSP_WAVELET_DN removes the low frequency part of the
%   signal and solves a l1 minimization problem to remove the noise.
%
%   Figure 1: Result of filtering
%
%      We observe that the wavelet denoising allows rapid change on the
%      graph signal. This is not possible with a simple low pass filtering.
%
%   Figure 2: Choosen wavelet filterbank
%
%      Here we use a mexican hat frame for the wavelet construction. 
%   
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/demos/gsp_demo_wavelet_dn.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

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

