function gsp_plot_signal_spectral(G,f_hat,param)
%GSP_PLOT_SIGNAL_SPECTRAL  Plot a graph signal in the graph spectral domain
%   Usage:  gsp_plot_signal_spectral(G,signal);
%           gsp_plot_signal_spectral(G,signal,param);
%
%   Input parameters:
%         G          : Graph or laplacian eigenvalues
%         f_hat      : Graph signal in the spectral domain.
%         param      : Optional variable containing additional parameters.
%   Output parameters:
%         none
%
%   'gsp_plot_signal_spectral(G,signal)' plots a graph signal in the graph
%   Fourier domain. 
%
%   Warning: to use this function, the Fourier basis of the graph should be
%   first computed. You can do it with::
%
%       G = gsp_full_eigen(G);
%
%   Example:::
%
%       N = 32;
%       G = gsp_path(N);
%       G = gsp_compute_fourier_basis(G);
%       f = sin((1:N)'*2*pi*4/N);
%       fhat = gsp_gft(G,f);
%       gsp_plot_signal_spectral(G,fhat);
%
%   Additional parameters
%   ---------------------
%
%   * *param.highlight* : Highlight one frequency component.
%   * *param.plot_abs*  : Option to plot the absolute value of f_hat
%
%   See also:  gsp_plot_graph gsp_plot_signal
%

% Author : David I Shuman, Nathanael Perraudin
% Testing: test_plotting

  
% Read input parameters
if nargin < 3
    param = struct;
end

if ~isfield(param, 'highlight'), param.highlight = 0; end
if ~isfield(param, 'plot_abs'), param.plot_abs = 1; end

if param.plot_abs
    f_hat=abs(f_hat);
end

if isstruct(G)
   E = G.e; 
end

cla;
hold on;
stem(E,f_hat,'LineWidth',1)
if param.highlight>0
    stem(E(highlight),f_hat(highlight),'filled','r','LineWidth',1);
end
tic_int=ceil(max(E)/10);
set(gca,'XTick',0:tic_int:ceil(max(E)));
set(gca,'color','none')

if sign(min(f_hat))>0
    aymin = 0.9*min(f_hat);
else
    aymin = 1.1*min(f_hat);
end
aymax = max(1.1*max(f_hat),1.1*min(f_hat)+eps);

axis([-.1 ceil(max(E)) aymin aymax]);
% set(gca,'fontsize',18)
box on;
hold off;
end
