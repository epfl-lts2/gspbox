%GSP_DEMO_PYRAMID Demonstration of the use of the graph pyramid
%
%   In this demonstration file, we show how to reduce a graph using the
%   GSPBox. Then we apply the pyramid to a simple signal.
%
%   The function GSP_GRAPH_MULTIRESOLUTION computes the graph pyramid for you:
%
%             param.sparsify = 1;
%             Gs = gsp_graph_multiresolution(G, Nlevel,param);
%
%   Gs is a cell array of graph representing the pyramid of graphs. Here
%   all optional parameter are important:
%
%    param.sparsify*: When a graph is reduced, the density of edges has
%     tendency to be amplified. One way to counterbalance this effect is to
%     sparsify the graph for each sublevel. The function
%     GSP_GRAPH_SPARSIFY is used to perform this operation. However, this
%     could lead to bad graphs (disconnected for instance).
%    sparsify_epsilon*: Parameter epsilon used in the sparsification
%    param.filters*: is a cell array of filters (or a single filter).
%     Thoses filter will be used in the analysis and synthesis operator.
%
%   Let's display the results:
%
%             figure;
%             for ii = 1:numel(Gs)
%                 subplot(2,3,ii)
%                 gsp_plot_graph(Gs{ii})
%                 title(['Reduction level: ', num2str(ii-1)]);
%             end
%
%   Figure 1: Reduction of the graph
%
%      
%
%   Now that we have precomputed the pyramid of graphs, we can apply it to
%   a signal. Here we create a signal that is smooth over the graph but
%   with a big discontinuity in the middle
%
%   Figure 2: Signal to be analysed
%
%
%
%   The graph pyramid can be simply applied thanks to the function 
%   GSP_PYRAMID_ANALYSIS:
%
%             [ca,pe]=gsp_pyramid_analysis(Gs,f);
%
%   ca contains the coarse approximation of each level and pe the
%   prediction errors. Let's display them:
%
%             figure
%             paramplot.show_edges = 0;
%             for ii = 1:numel(Gs)
%                 subplot(2,3,ii)
%                 gsp_plot_signal(Gs{ii},pe{ii},paramplot);
%                 title(['P. E. level: ', num2str(ii-1)]);
%             end
% 
%             figure
%             for ii = 1:numel(Gs)
%                 subplot(2,3,ii)
%                 gsp_plot_signal(Gs{ii},ca{ii},paramplot)
%                 title(['C. A. level: ', num2str(ii-1)]);
%             end
%
%   Figure 3: Prediction errors
%
%
%
%   Figure 4: Coarse approximations
%
%
%   Finally, you can perform a synthesis operation using the function
%   GSP_PYRAMID_SYNTHESIS :
%   
%             coeff = gsp_pyramid_cell2coeff(ca,pe);
%             f_pred = gsp_pyramid_synthesis(Gs,coeff);
%
%   The function GSP_PYRAMID_CELL2COEFF remove all unecessary
%   coefficients and keep only the last sublevel and the prediction error.
%
%   Enjoy!
%
%   Url: http://lts2research.epfl.ch/gsp/doc/demos/gsp_demo_pyramid.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

clear;
close all;

%% Parameters


% Load a graph
N=512;
param_sensor.distribute = 1;
G = gsp_sensor(N,param_sensor);

% Number of level
Nlevel = 5;

%% Compute the graph decomposition

param.sparsify = 1;
Gs = gsp_graph_multiresolution(G, Nlevel,param);
Gs = gsp_compute_fourier_basis(Gs);
% % Compute lmax for each subgraph
%  Gs = gsp_estimate_lmax(Gs);

% Display the results
figure;
for ii = 1:numel(Gs)
    subplot(2,3,ii)
    gsp_plot_graph(Gs{ii})
    title(['Reduction level: ', num2str(ii-1)]);
end

%% Apply the pyramid to a signal

% create a signal
f = ones(N,1);
f(1:N/2) = -1;
f = f+10*Gs{1}.U(:,8);

figure
gsp_plot_signal(G,f)


[ca,pe]=gsp_pyramid_analysis(Gs,f,Nlevel);

figure
paramplot.show_edges = 1;
for ii = 1:Nlevel
    subplot(2,3,ii)
    gsp_plot_signal(Gs{ii},pe{ii},paramplot);
    title(['P. E. level: ', num2str(ii-1)]);
end

figure
for ii = 1:Nlevel+1
    subplot(2,3,ii)
    gsp_plot_signal(Gs{ii},ca{ii},paramplot)
    title(['C. A. level: ', num2str(ii-1)]);
end



%% Perform synthesis

f_pred = gsp_pyramid_synthesis(Gs,ca{end},pe);

err = norm(f_pred-f)/norm(f);

fprintf('The relative reconstruction error is : %g \n',err);

