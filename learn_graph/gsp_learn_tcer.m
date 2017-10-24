function [TCER, TCER_rand] = gsp_learn_tcer(G, X, params)
%GSP_LEARN_TCER Learn the Total Cumulative Energy Residual
%   Usage: TCER = gsp_learn_tcer(G, X);
%          TCER = gsp_learn_tcer(G, X, params);
%          [TCER, TCER_rand] = gsp_learn_tcer(...);
% 
%   Input parameters:
%       G           : the graph
%       X           : a data matrix
%   Output parameters: 
%       TCER        : the computed total commulative energy residual for G 
%       TCER_rand   : the computed total commulative energy residual for a random basis 
% 
%   Total Cumulative Energy Residual (TCER): This is a number from [0 to 1]
%   describing how well a graph fits a given data matrix X, or
%   distribution.
%
%   Optional parameters
%   -------------------
%
%   * *params.s*      : s = svd(X);
%   * *params.verbose*: 0 = nothing, 1 = plot cum energy, 2 = plot against
%     random basis as a baseline (default 0)
%   * *params.sort*   : sort basis columns to get the minimum residual
%     possible (default 0)
% 
%   See also: gsp_good_graph_index, gsp_stationarity_ratio
% 


% Author  : Andreas Loukas, Vassilis Kalofolias
% Date    : 15 Nov 2016

% Handle input
if nargin < 3, params = struct(); end
if not(isfield(params, 'verbose')); params.verbose = 0; end;
if not(isfield(params, 'sort'));    params.sort = 0;    end;
if not(isfield(params, 's'));       params.s = svd(full(X));  end;

s = params.s;

if numel(s)< G.N
    s = [s; zeros(G.N-numel(s),1)];
end

if not(isfield(G, 'U'))
    G = gsp_compute_fourier_basis(G);
end

% frobenius norm of X
norm_X_fro = norm(s);

% sorting already done by svd.
% var_cumsum_svd = cumsum(sort(s.^2, 'descend')) ./ norm_X_fro^2;
var_cumsum_svd = cumsum(s.^2) ./ norm_X_fro^2;

% frobenius norm of X
norm_X_fro = norm(X, 'fro');

% for the graph
if params.sort
    var_cumsum_s_G = cumsum(sort(sum(abs(G.U' * X).^2, 2) ./ norm_X_fro^2, 'descend'));
else
    var_cumsum_s_G = cumsum(sum(abs(G.U' * X).^2, 2) ./ norm_X_fro^2);
end

% area under the curve
AUC_SVD   = sum(var_cumsum_svd);
AUC_G     = sum(var_cumsum_s_G);
TCER      = 1 - AUC_G/AUC_SVD;

% compute also the result for a random basis 
AUC_RB    = sum(linspace(1/G.N, 1, G.N));
TCER_rand = (AUC_SVD - AUC_RB)/(AUC_SVD);

if params.verbose >= 1
    
    figure; hold on;

    plot(var_cumsum_svd)
    plot(var_cumsum_s_G)
    [Urand, ~] = qr(randn(G.N));        % random basis
    var_cumsum_s_rand = cumsum(sum(abs(Urand' * X).^2, 2) ./ norm_X_fro^2);
    plot(var_cumsum_s_rand)

    xlabel('k (number of basis vectors used)');
    ylabel('cumulative variance');
    legend('svd', 'G', 'random basis')    
end

% plot good quality residual figure
if params.verbose >= 2
    figure; 
    plot(var_cumsum_svd);
    hold on;
    %     plot(cum_energy_graph);
    plot(var_cumsum_s_G, 'color', [0.8500    0.3250    0.0980]);
    n_nodes = length(var_cumsum_svd);
    xf = [1:n_nodes, fliplr(1:n_nodes)];
    yf = [zeros(n_nodes,1); flipud(var_cumsum_svd)];
    fill(xf, yf, [     0    0.4470    0.7410],'facealpha', .2);

    xf = [1:n_nodes, fliplr(1:n_nodes)];
    yf = [zeros(n_nodes,1); flipud(var_cumsum_s_G)];
    fill(xf, yf, [0.8500    0.3250    0.0980],'facealpha', .2);

    plot(var_cumsum_s_G, 'color', [0.8500    0.3250    0.0980]);
    plot(var_cumsum_svd, 'color', [     0    0.4470    0.7410]);
    legend('Ground truth covariance', 'Graph', 'Ground truth total C.E.', 'Graph total C.E.','location', 'se')
%     title(['graph from ', num2str(k), ' samples'])
    ylabel('Expected cumulative energy');
    text(n_nodes/2, .5, ['Residual: ', num2str(100*TCER, 3), '%']);
    xlabel('k');
    ylim([0,1]);
end