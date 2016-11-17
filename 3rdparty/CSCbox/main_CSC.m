% This is the main function where you enter parameters, your graph to analyze, 
% and launch the CSC algorithm
%
% Copyright (C) 2016 Nicolas Tremblay, Gilles Puy.
% This file is part of the CSCbox (Compressive Spectral Clustering toolbox)
%
% The CSCbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The CSCbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%     N. Tremblay, G. Puy, R. Gribonval and P. Vandergheynst.
%     Compressive Spectral Clustering.
%     ArXiv e-prints:1602.02018 Feb. 2016.

%% ====================== PARAMETERS ===========================
%%%

param_CSC.lap_type = 'normalized'; % type of Laplacian used; 'normalized' or 'combinatorial'
param_CSC.n_factor = 2; % n = param_CSC.n_factor * G.k * log(G.k);
param_CSC.d_factor = 4; % d = param_CSC.d_factor * log(n);
param_CSC.sampling = 'uniform'; % this is the sampling distrib of nodes. One may also try 'weighted'. In the paper, the sampling distribution is uniform. 
param_CSC.poly_order = 50; % polynomial approximation order, p in the paper
param_CSC.regu = 1e-3; % regularisation parameter of the solver, \gamma in the paper
param_CSC.solver = 'gmres'; % the solver; either 'gmres' or 'cgs'

%%%
%% ====================== GRAPH INPUT ===========================
%%%
% at the end of this "GRAPH INPUT" code section, one needs 
% - G.k, the number of classes to find; 
% - G.W, the adjacency matrix of the graph (should be sparse and symmetrical); 
% - G.N, the number of nodes in the graph; 
% - and G.truth, if it exists, the ground truth of the clustering.

%% ------------------- a possible choice: the SBM benchmark
G.N = 1e4; % number of nodes
G.k = 100;  % number of communities
G.c = 16; % average degree

%---- if homogeneous com size
com_size=round(ones(1,G.k).*(G.N/G.k));

%---- if heterogeneous com size (eg with N=10^3)
% com_size=[5,10,15,20,25,30,35,40,45,50,50,55,60,65,70,75,80,85,90,95];
% epsi_c=0.20;

%---- difficulty of community detection
epsi_c=(G.c-sqrt(G.c))/(G.c+sqrt(G.c)*(G.k-1)); %maximum difficulty
epsi=epsi_c/4; % the closer to espi_c, the more difficult is the clust. task

%---- create SBM graph
G = create_SBM(G.N,G.k,G.c,epsi,com_size);

% % %% ------------------- another possible choice: the amazon graph
% % 
% %---- Load data
% edges = importdata('com-amazon.ungraph2.txt');
% %---- Get set of identifiers appearing in the above file
% nodes = unique(edges(:));
% %---- Map nodes to 1:N and recompute edges accordingly
% [~, edges(:, 1)] = ismember(edges(:, 1), nodes);
% [~, edges(:, 2)] = ismember(edges(:, 2), nodes);
% %---- Build adjacency matrix
% G.N = numel(nodes);
% G.W = sparse(edges(:, 1), edges(:, 2), 1, G.N, G.N);
% %---- Ensure matrix is symmetric
% G.W = G.W + G.W'; G.W(G.W>0) = 1;
% clear edges nodes;
% %---- Nb. Clusters
% G.k = 250;

%%%
%% ====================== perform CSC ===========================
%%%

[C_est, lk_est, time_CSC, IDX_LD, ind_obs, weight_VD] = CSC(G,param_CSC);
% C_est is the G.N x G.k interpolated cluster indicator functions
% lk_est is the estimated k-th eigenvalue of the Laplacian
% time_CSC records all times of computation
% IDX_LD is the low-dimensional k-means result
% ind_obs are the sampled node numbers
% weight_VD is the estimated variable density (it was actually used only if
% param_CSC.sampling is equal to 'weighted'). 

C_est=C_est./repmat(sqrt(sum(C_est.^2)),G.N,1); %normalize
[~, IDX_CSC] = max(C_est, [], 2);
% that's because we look for hard clustering here. If one looks for fuzzy
% clustering, C_est contains more information than IDX_CSC

% compute performance if the ground truth exists:
perf_CSC=PartAgreeCoef_ARonly(G.truth,IDX_CSC);
fprintf('\nCSC: Performance (ARI) = %6.3f\n', perf_CSC); % ARI= Adjusted Rand Index
fprintf('CSC: Time of computation (s) = %6.3f\n', time_CSC.total);

% % or else compute modularity to have some "community measure" of the result:
% modu_CSC=compute_modularity(IDX_CSC,G.W);
         
%%%
%% ================== perform SC (Ng et al.)====================
%%%

[IDX_SC, G.Uk, G.Dk, time_SC] = SC(G,param_CSC.lap_type);
% IDX_SC is the clustering result of classical SC
% G.Uk contains the k first eigenvectors of the Laplacian
% G.Dk contains the k first eigenvalues
% time_SC records the time of comp. of the different steps of the algorithm

lk_true=G.Dk(end); %only to record the "true" value of lambda_k

% compute performance if the ground truth exists:
perf_SC=PartAgreeCoef_ARonly(G.truth,IDX_SC);
fprintf('\nSC: Performance (ARI) = %6.3f\n', perf_SC); % ARI= Adjusted Rand Index
fprintf('SC: Time of computation (s) = %6.3f\n', time_SC.total);

% % or else compute modularity to have some "community measure" of the result:
% modu_SC=compute_modularity(IDX_SC,G.W);

