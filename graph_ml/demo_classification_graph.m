
%% Initialization
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/demo_classification_graph.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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

clear 
close all;

%% Handle the data
% USPS dataset
[X_obs, labels_obs, X_unobs, labels_unobs] = load_usps_full();
X_obs = X_obs';
X_unobs = X_unobs';

N_obs = 1000;
N_unobs = 1000;

ind_kept_obs = randperm(size(X_obs, 1), N_obs);
ind_kept_unobs = randperm(size(X_unobs, 1), N_unobs);

X_obs = X_obs(ind_kept_obs, :);
labels_obs = labels_obs(ind_kept_obs);
X_unobs = X_unobs(ind_kept_unobs, :);
labels_unobs = labels_unobs(ind_kept_unobs);

[labels_obs, ind] = sort(labels_obs);
X_obs = X_obs(ind, :);

[labels_unobs, ind] = sort(labels_unobs);
X_unobs = X_unobs(ind, :);

% sizes
N_obs = size(X_obs, 1);
N_unobs = size(X_unobs, 1);

% % Synthetic dataset
% % sizes
% Nl = 10;
% Nnl = 100;
% x = rand(Nl,2);
% xx = rand(Nnl,2);
% y = double(x(:,1)<0.5);
% yy = double(xx(:,1)<0.5);



%% Define the error function

err = @(x) sum(sum(abs(x((N_obs+1):end)-labels_unobs)>0))/N_unobs;


%% Create the graphs
param.use_flann = 1;
param.k = 6;

G = gsp_nn_graph([X_obs;X_unobs],param);
G = gsp_create_laplacian(G,'normalized');


% Make the assotiated mask
M = zeros(N_obs+N_unobs,1);
M(1:N_obs) = 1;

% Starting point
y_start = [labels_obs; zeros(N_unobs,1)]';

%% Method 1: TIK 

s_tik = gsp_classification_tik(G,M,y_start, 0 , param);
err_tik = err(s_tik);

%% Method 2:  TV 
param.maxit = 2000;
s_tv = gsp_classification_tv(G,M,y_start, 0 , param);
err_tv = err(s_tv);

%% Method 3:  TV 
s_tv_new = gsp_classification_tv_new(G,M,y_start, 0 , param);
err_tv_new = err(s_tv_new);

%% Plotting: only for synthetic dataset
% 
% figure(1)
% 
% subplot(221);
% gsp_plot_signal(G,[labels_obs;labels_unobs]);
% title('Original signal')
% subplot(222);
% gsp_plot_signal(G,M);
% title('Mask')
% 
% subplot(223);
% gsp_plot_signal(G,s_tik);
% title('Tikonov classifier')
% 
% subplot(224);
% gsp_plot_signal(G,s_tv);
% title('TV classifier')


%%
figure; plot(s_tv(N_obs+1: end)); hold on; plot(labels_unobs)
figure; plot(s_tik(N_obs+1: end)); hold on; plot(labels_unobs)
figure; plot(s_tv_new(N_obs+1: end)); hold on; plot(labels_unobs)


