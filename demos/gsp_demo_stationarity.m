%GSP_DEMO_STATIONARITY Demonstration to use stationarity
%
%   This demonstration uses the theory developed in the paper "Stationary
%   graph signal processing?, Nathanael Perraudin and Pierre Vandergheynst,
%   arXiv:1601.02522, http://arxiv.org/pdf/1601.02522
%
%   We use Wiener optimization to solve an in-painting problem. We suppose
%   that the Power Spectral Density (PSD) of the input signal is unknown
%   and we estimate it using a single signal. The PSD is estimated using
%   the function |gsp_estimate_psd| ::
%
%             parampsd.Nfilt = 100;
%             psd_est = gsp_estimate_psd(G,X1,parampsd);
%
%   The parameter *param.Nfilt* selects the number of windows to use for
%   the PSD estimation. A large number leads to a finer approximation.
%   However it 1) increases the computation time for the estimation 2)
%   necessitates a larger order of Chebyshev polynomial if exact filtering
%   is not used. the Chebyshev order can be changed with::
% 
%             parampsd.order = 50; 
%
%   Once the PSD is computed, one can estimate a signal with the function
%   |gsp_wiener_inpainting| ::
%
%             S = gsp_wiener_inpainting(G, Y, Mask, psd_est, sigma^2,param);
%
%   This function solves the Wiener optimization problem thanks to the
%   UNLocBoX.  You might want to change the precision of the recovery. To do
%   so, you can set the parameter *maxit* and *tol* of the structure
%   *param*::
%
%             param.maxit = 1000; % Maximum number of iterations
%             param.tol = 1e-8; % Tolerance to stop iterating
%
%   Figures 1,2 and 3 show an inpainting example. Figure 4 shows the
%   estimation of the PSD from only one measurement.
%
%   .. figure::
%
%      We generate a stationary signal
%
%       
%
%   .. figure::
%
%      We randomly select 50% of the signal
%
%
%
%   .. figure::
%
%      We reconstruction the signal using Wiener optimization
%
%       
%
%   .. figure::
%
%      True VS approximated PSD 
%
%       
%
%
%   References: perraudin2016stationary
 
 
 
% Author : Nathanael Perraudin
% Date: 25 Mai 2016
 
%% Initialization
clear
close all;
gsp_reset_seed(0)
 
% Parameters
N = 400; % Size of the graph
K = 1; % Number of signal(s) to estimate the PSD
g = @(x) sin(3*x).*exp(-3*x); % generative filter
psd = @(x) (g(x)).^2; % psd
sigma = 0; % noise parameter
 
 
 
%% Generate data and create the problem
 
% Create a graph
G = gsp_sensor(N);
% Compute the Fourier basis (if the graph is small)
G = gsp_compute_fourier_basis(G);
% G = gsp_estimate_lmax(G); % if the graph is large
 
% Create the data
w = randn(N,K);
X1 = gsp_filter_analysis(G,g,w);
w2 = randn(N,1);
X2 = gsp_filter_analysis(G,g,w2);
 
% Create a mask
Mask = double(rand(N,1)>0.5); % Mask
% Create the measurements
Y = Mask.*X2 + sigma *randn(N,1);
Y(logical(1-Mask)) = min(Y(:));
 
%% Solve the problem in two steps
 
% 1) Estimate the PSD
parampsd.Nfilt = 100;
% parampsd.order = 50; % to increase the order of the Chebyshev polynomial
psd_est = gsp_estimate_psd(G,X1,parampsd);
 
% 2) Prediction
S = gsp_wiener_inpainting(G, Y, Mask, psd_est, sigma^2);
 
% % Compute the error
% snr_pred = snr(X2(logical(1-Mask)),S(logical(1-Mask)));
 
%% Display the results
 
figure(1); gsp_plot_signal(G,X2); title('Original signal');
figure(2); gsp_plot_signal(G,Y); title('Measurments');
figure(3); gsp_plot_signal(G,S); title('Wiener inpainting');
 
figure(4); 
paramplot.plot_eigenvalues = 0;
paramplot.cla = 0;
gsp_plot_filter(G,psd,paramplot);
hold on
gsp_plot_filter(G,psd_est,paramplot);
legend('True PSD','Estimate PSD')