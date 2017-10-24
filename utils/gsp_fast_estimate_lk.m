function [lambda_k, info] = gsp_fast_estimate_lk(G, k, param)
%GSP_FAST_ESTIMATE_LK Estimation of k-th eigenvalue of any graph Laplacian
%   Usage:  lambda_k = gsp_fast_estimate_lk(G,k);
%           [lambda_k, info] = gsp_fast_estimate_lk(G, k, param);
%
%   Input parameters :
%         G          : Graph structure.
%         k          : Index of the eigenvalue to estimate.
%         param      : Optional parameters
%   Output parameters:
%         lambda_k   : Value of the $k$ th eigenvalue
%         info       : Statistics for each estimation ('lk_est', 'k_est', 'calls')
%
%   'gsp_fast_estimate_lk(G,k)' computes an estimation of the $k$ th
%   eigenvalue of the Laplacian of G using accelerated eigencount technique
%   assuming local uniformity over the spectrum as described in the
%   reference below.
%
%   Example:::
%
%      G = gsp_sensor(256);
%      G = gsp_estimate_lmax(G);
%      k = 8;
%      param.order = 100;
%      lambda_k = gsp_fast_estimate_lk(G, k, param);
%       
%
%   Additional parameters
%   ---------------------
%  
%   * *param.filter*  : Select the filter to be used for the computation. 
%     * 'lp-ch'   : Chebyshev polynomial approximation
%     * 'lp-jch'  : Jackson-Chebyshev polynomial approximation
%     * 'expwin'  : Exponentially decreasing polynomial approximation. Default: 'lp-jch'
%   * *param.order* : Degree of the polynomial approximation (default=50).
%   * *param.nb_estimation* : Number of estimations to average.
%   * *param.nb_features* : Number of features to filter.
%   * *param.max_calls* : Max number of calls for an estimation.
%   * *param.hint_lambda_max* : Hint on upper bound of lk.
%   * *param.hint_c_max* : Number of eigenvalues up to hint_lmax.
%   * *param.hint_lambda_min* : Hint on lower bound of lk.
%   * *param.hint_c_min* : Number of eigenvalues up to hint_lmin.
%   * *param.verbose* : Verbosity level (0 no log - 1 display warnings) (default 1).
%    
% 
%   References: paratte2016fast
%

% Author: Johan Paratte, Lionel Martin
% Date: 3 November 2016


if nargin < 3, param = struct; end
if ~isfield(param, 'nb_estimation'), param.nb_estimation = 1; end
if ~isfield(param, 'nb_features'), param.nb_features = 2*round(log(G.N)); end
if ~isfield(param, 'max_calls'), param.max_calls = 10; end
if ~isfield(param, 'filter'), param.filter = 'lp-jch'; end
if ~isfield(param, 'order'), param.order = 50; end
if ~isfield(param, 'verbose'), param.verbose = 1; end

if ~isfield(G, 'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_FAST_ESTIMATE_LK: The variable lmax is not available.', ...
            'The function will compute it for you. However, if you apply ', ...
            'many time this function, you should precompute it using the ', ...
            'function: gsp_estimate_lmax.']);
end

if ~isfield(param, 'hint_lambda_max') || ~isfield(param, 'hint_c_max')
    param.hint_lambda_max = G.lmax;
    param.hint_c_max = G.N;
end
if ~isfield(param, 'hint_lambda_min') || ~isfield(param, 'hint_c_min')
    param.hint_lambda_min = 0;
    param.hint_c_min = 0;
end


% List of estimations for lambda_k
lambda_k_est = zeros(param.nb_estimation, 1);
k_est = zeros(param.nb_estimation, 1);
calls = zeros(param.nb_estimation, 1);

% Perform nb_estimation on different sets of feature vectors
for ind_est = 1:param.nb_estimation
    
    % Random signals (fixed for one estimation)
    Sig = randn(G.N, param.nb_features)/sqrt(param.nb_features);
    
    % Initial values
    lmin = param.hint_lambda_min;
    lmax = param.hint_lambda_max;
    cmin = param.hint_c_min; % eigencount above lmin
    cmax = param.hint_c_max; % eigencount above lmax

    l_est = lmin + (k-cmin)*(lmax - lmin)/(cmax-cmin);
    [lambda_k_est(ind_est), k_est(ind_est), calls(ind_est)] = rec_estimate(G, k, Sig, l_est, lmin, lmax, cmin, cmax, param);
end

% Final estimation
lambda_k = mean(lambda_k_est);
info = struct('lk_est', lambda_k_est, 'k_est', k_est, 'calls', calls);

end

function [lk, count, nb_calls] = rec_estimate(G, k, signal, lest, lmin, lmax, cmin, cmax, lbest, cbest, nb_calls, param)
    if nargin == 9
        param = lbest;
        nb_calls = 0;
        if (abs(cmin - k) < abs(cmax - k))
            lbest = lmin;
            cbest = cmin;
        else
            lbest = lmax;
            cbest = cmax;
        end
    end
    
    nb_calls = nb_calls + 1;

    switch param.filter
        case 'lp-ch'
            [pcoefs, ~] = jackson_cheby_poly_coefficients(0, lest, [0, G.lmax], param.order);

        case 'lp-jch'
            [~, pcoefs] = jackson_cheby_poly_coefficients(0, lest, [0, G.lmax], param.order);

        case 'expwin'
            ew = gsp_design_expwin(G, lest/G.lmax);
            pcoefs = gsp_cheby_coeff(G, ew, param.order);

        otherwise
            disp('Unknown filter type');     
    end

    X = gsp_cheby_op(G, pcoefs, signal);
    c_abs = round(sum(X(:).^2));

    ck = k - cmin; % eigenvalues missing above lmin
    
    cl = c_abs - cmin; % eigenvalues found in [lmin, lest]
    cr = cmax - c_abs; % eigenvalues remaining in [lest, lmax]
    
    if abs(k - c_abs) <= abs(k - cbest)
        cbest = c_abs;
        lbest = lest;
    end
    
    if param.verbose
        fprintf('Iter #%d - lest: %f [%d] - lmin: %f [%d] - lmax: %f [%d] || best: %f [%d]', nb_calls, lest, c_abs, lmin, cmin, lmax, cmax, lbest, cbest);
        if (cl == 0 || cr == 0), fprintf('   [DICHOTOMY]\n'); else fprintf('\n'); end
    end
    
    if cbest == k
       lk = lbest;
       count = cbest;
       return;
    end
    
    %Special cases : c_abs = c_min or c_abs = c_max
    if cl == 0 || cr == 0
        if cl == 0
            lmin = lest;
        elseif cr == 0
            lmax = lest;
        end

        lest = (lmax + lmin) / 2;
    else
        if cl > ck
            lmax = lest;
            lest = lmin + ck * (lest - lmin)/cl;
            cmax = c_abs;
        end

        if cl < ck
            lmin = lest;
            lest = lest + (ck - cl)*(lmax - lest)/cr;
            cmin = c_abs;
        end
    end
    
    if nb_calls == param.max_calls
       lk = lbest;
       count = cbest;
       return;
    end

    [lk, count, nb_calls] = rec_estimate(G, k, signal, lest, lmin, lmax, cmin, cmax, lbest, cbest, nb_calls, param);
end
