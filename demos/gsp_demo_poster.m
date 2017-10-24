%GSP_DEMO_POSTER A 6 steps demonstration
%    
%   This file contains a 6 steps demonstration of the GSPBOX
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/demos/gsp_demo_poster.html

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
