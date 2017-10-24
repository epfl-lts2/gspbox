function [  ] = gsp_plot_sgram( G,A,param )
%GSP_PLOT_SGRAM Plot graph spectrogram
%   Usage:  gsp_plot_sgram( G,A );
%           gsp_plot_sgram( G,A,param );
%
%   Input parameters:
%         G     : Graph
%         A     : Graph windowed Fourrier transform
%         param : Structure of optional parameter
%   Output parameters:
%         none
%
%   *param* is a structure of optional parameter with
%
%   * *param.colorbar*: Use the colorbar (default 1)
%
%   Example:::
%
%           N = 15;
%           G = gsp_ring(2*N);
%           G = gsp_compute_fourier_basis(G);
%           x = [0:N,(N-1):-1:1]';
%           s = 3;
%           g = exp(-(x-1).^2/s^2);
%           f = gsp_localize(G,g,N);
%           c = gsp_gwft(G,f,g);
%           gsp_plot_sgram(G,c);
%   
%   See also: gsp_plot_signal, gsp_plot_graph, gsp_plot_signal_spectral
%


% Author: Nathanael Perraudin
% Date  : 09.12.2013
% testing: test_plotting

% Optional parameter handling
if nargin<3
    param=struct;
end

if ~isfield(param, 'colorbar'), param.colorbar = 1; end;

imagesc(1:size(A,2), 0:size(A,1)-1,abs( A));

% Hack to overpass a matlab bug with latex interpretex
latex = get(gca,'DefaultTextInterpreter');
set(gca,'DefaultTextInterpreter','Tex');

xlabel('Nodes');
ylabel('Freqencies');

set(gca,'DefaultTextInterpreter',latex);

if param.colorbar
    colorbar
end




end

