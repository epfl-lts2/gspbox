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
%   param is a structure of optional parameter with
%
%    param.colorbar*: Use the colorbar (default 1)
%
%   Example:
%
%           N = 15;
%           G = gsp_ring(2*N);
%           G = gsp_compute_fourier_basis(G);
%           x = [0:N,(N-1):-1:1]';
%           s = 3;
%           g = exp(-x.^2/s^2);
%           f = gsp_modulate(G,gsp_translate(G,g,N),N);
%           c = gsp_gwft(G,f,g);
%           gsp_plot_sgram(G,c);
%   
%   See also: gsp_plot_signal, gsp_plot_graph, gsp_plot_signal_spectral
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/plotting/gsp_plot_sgram.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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


