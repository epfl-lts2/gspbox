% function [G, lambda_k, cum_coh, X] = estimate_lambda_k(G, k, param)
% 
% This function estimates the k-th eigenvalue of the Laplacian matrix of a graph G.
% Inputs:
% - G.L should be the sparse symmetrical Laplacian matrix of the graph.
% - G.lmax should be the maximal eigenvalue of G.L.
% - k should be the index of the eigenvalue you want to estimate. 
% - param are parameters. 
% *param.nb_estimation is the number of estimation. Default is 1. 
% *param.nb_features is the number of randomp features for estimation. Default is 2*round(log(G.N)). 
% *param.epsilon is the precision one desires. Default is 1e-1. 
% *param.hint_lambda_max is a hint given to the algorithm to start the dichotomy 
% with param.hint_lambda_max as upper bound instead of G.lmax. Default is G.lmax. 
% *param.hint_lambda_min is a hint to start the dichotomy with that lower bound 
% instead of 0. Default is 0. 
% *param.jackson=0 if one wants to use Chebychev polynomials. param.jackson=1 if 
% one wants to use Jackson-Chebychev polynomials. Default is 1.
% *param.order is the order of the polynomial filters used. Default is 50.
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

function [G, lambda_k, cum_coh, X] = estimate_lambda_k(G, k, param)


% Parameters
if nargin < 3
    param = struct;
end
if ~isfield(param, 'nb_estimation'),
    param.nb_estimation = 1;
end
if ~isfield(param, 'nb_features'),
    param.nb_features = 2*round(log(G.N)); %k;
end
if ~isfield(param, 'epsilon'),
    param.epsilon = 1e-1;
end
if ~isfield(param, 'hint_lambda_max'),
    param.hint_lambda_max = G.lmax;
end
if ~isfield(param, 'hint_lambda_min'),
    param.hint_lambda_min = 0;
end
if ~isfield(param, 'jackson'),
    param.jackson = 1;
end
if ~isfield(param, 'order'),
    param.order = 50;
end

% List of estimations for lambda_k
norm_Uk = zeros(G.N, param.nb_estimation);
lambda_k_est = zeros(param.nb_estimation, 1);

% Perform nb_estimation on different of set feature vectors
for ind_est = 1:param.nb_estimation
    
    % Random signals (fixed for one estimation)
    Sig = randn(G.N, param.nb_features)*1/sqrt(param.nb_features);
    
    % Search by dichotomy
    counts = 0;
    lambda_min = param.hint_lambda_min;
    lambda_max = param.hint_lambda_max;
    while (counts~=k) && ((lambda_max - lambda_min)/lambda_max > param.epsilon)
        % Middle of the interval
        lambda_mid = (lambda_min+lambda_max)/2;
        % Filter
        [ch, jch] = jackson_cheby_poly_coefficients(0, lambda_mid, [0 G.lmax], param.order);
        % Filtering
        if param.jackson==2
            [X,Xjch] = gsp_cheby_op2(G, ch(:), jch(:), Sig);
            countsch = round(sum(X(:).*Sig(:)));
            countsjch = round(sum(Xjch(:).*Sig(:)));
            counts=round((countsch+countsjch)/2);
        elseif param.jackson==1
            X = gsp_cheby_op(G, jch(:), Sig);
            counts = round(sum(X(:).^2));
            %counts = round(sum(X(:).*Sig(:)));
            %counts=round((sum(X(:).^2)+sum(X(:).*Sig(:)))/2);
        else
            X = gsp_cheby_op(G, ch(:), Sig);
            counts = round(sum(X(:).^2));
            %counts = round(sum(X(:).*Sig(:)));
            %counts=round((sum(X(:).^2)+sum(X(:).*Sig(:)))/2);
        end
        % Check results
        
        if counts>k, lambda_max = lambda_mid;
        elseif counts<k lambda_min = lambda_mid;
        end
    end
    
    % Store result
    lambda_k_est(ind_est) = (lambda_min+lambda_max)/2;
    
    norm_Uk(:, ind_est) = sum(X.^2, 2);
end

% Final estimation
G.lk = mean(lambda_k_est);
lambda_k = G.lk;
cum_coh = mean(norm_Uk, 2);

end
