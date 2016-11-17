function [ basis, filt_sig, param ] = gsp_eigenspace_estimation( G, k, param )
%GSP_EIGENSPACE_ESTIMATION Estimation of first eigenvectors of any graph Laplacian
%   Usage:  basis = gsp_eigenspace_estimation(G,k);
%           [basis, approx_U, param] = gsp_eigenspace_estimation(G,k,param);
%
%   Input parameters :
%         G          : Graph structure.
%         k          : Dimension of the subspace.
%         param      : Optional parameters
%   Output parameters:
%         basis      : Approximated basis of k first eigenvectors
%         filt_sig   : Filtered random signal
%         param      : Optional parameters (with new entries)
%
%   'gsp_eigenspace_estimation(G,k)' computes an estimation of the first 
%   k eigenvectors of the Laplacian of G using Gaussian random signal
%   filtering, following the FEARS method described in paratte2017fast.
%
%
%   Example:
%
%         G = gsp_sensor(256);
%         k = 8;
%         param.order = 100;
%         Uk_est = gsp_eigenspace_estimation(G, k, param);
%         G = gsp_compute_fourier_basis(G);
%         proj_energy = norm(Uk_est'  G.U(:, 1:k), 'fro');
%       
%
%   Additional parameters
%   ---------------------
%  
%    param.filter  : Select the filter to be used for the computation. 
%      'lp-ch'   : Chebyshev polynomial approximation
%      'lp-jch'  : Jackson-Chebyshev polynomial approximation
%      'expwin'  : Exponentially decreasing polynomial approximation Default: 'lp-jch'
%    param.order : Degree of the polynomial approximation (default=50).
%    param.lk_est_method : Select the version of lk estimation.
%      'fast'      : Accelerated method using local uniformity assumption
%      'std'  : Usual method using dichotomy all the time Default: 'fast'
%    param.R : Random matrix to use (of size N > d, d >= k)
%     (default: Gaussian(0, 1/k) of size Nxk)
%    param.pcoef : Polynomial coefficients if already known.
%    param.lk : Estimated value of lambda_k if already known.
%    param.verbose : Verbosity level (0 no log - 1 display warnings) (default 1).   
%
% 
%   References:
%     J. Paratte and L. Martin. Fast eigenspace approximation using random
%     signals. arXiv preprint arXiv:1611.00938, 2016.
%     
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/embedding/gsp_eigenspace_estimation.php

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

% Author: Johan Paratte, Lionel Martin
% Date: 3 November 2016

if nargin < 3, param = struct; end
if ~isfield(param, 'filter'), param.filter = 'lp-jch'; end
if ~isfield(param, 'order'), param.order = 50; end
if ~isfield(param, 'lk_est_method'), param.lk_est_method = 'fast'; end
if ~isfield(param, 'R'), param.R = randn(G.N, k)/sqrt(k); end
if ~isfield(param, 'verbose'), param.verbose = 1; end

assert(size(param.R, 1) == G.N && size(param.R, 2) >= k, 'The optional parameter R has wrong size.');

if ~isfield(param, 'pcoefs')
    if param.verbose, disp('Polynomial filtering required. Computing polynomial coefficients...'); end

    if ~isfield(G, 'lmax')
        G = gsp_estimate_lmax(G);
        warning(['GSP_EIGENSPACE_ESTIMATION: The variable lmax is not available.', ...
            'The function will compute it for you. However, if you apply ', ...
            'many time this function, you should precompute it using the ', ...
            'function: gsp_estimate_lmax.']);
    end

    if ~isfield(param, 'lk')
        if param.verbose, fprintf('Estimation of lambda_k');
        end

        tic;
        switch param.lk_est_method
            case 'fast'
                [param.lk, info] = gsp_fast_estimate_lk(G, k, param);
                if param.verbose, disp('using our accelerated method.'); end

            case 'std'
                [~, param.lk, ~, ~, nb_iter_lk, k_est_lk] = gsp_estimate_lk(G, k, param);
                if param.verbose, disp('using the standard method.'); end

            otherwise
                error('Unknown method for lk_est_method.');
        end
        t = toc;

        if param.verbose
            fprintf(['* Estimated lk: %d\n', ...
            '* Time to estimate lk: %f sec\n', ...
            '* in %d iterations with k_est=%d (target=%d)\n'], ...
            param.lk, t, mean(info.calls), mean(info.k_est), k);
        end
    else
        warning('lambda_k was provided to the method from param.');
    end

    tic;
    switch param.filter
        case 'lp-ch'
            [param.pcoefs, ~] = jackson_cheby_poly_coefficients(0, param.lk, [0, G.lmax], param.order);

        case 'lp-jch'
            [~, param.pcoefs] = jackson_cheby_poly_coefficients(0, param.lk, [0, G.lmax], param.order);

        case 'expwin'
            ew = gsp_design_expwin(G, param.lk/G.lmax);
            param.pcoefs = gsp_cheby_coeff(G, ew, param.order);

        otherwise
            error('Unknown filter type!');     
    end

    t = toc;
    if param.verbose, fprintf('* Time to compute polynomial coefficients: %f sec.\n', t); end

else
    if param.verbose, warning('pcoef was provided to the method from param.'); end
end

if param.verbose, disp('Filtering random signals...'); end
tic;
filt_sig = gsp_cheby_op(G, param.pcoefs, param.R);
t = toc;
if param.verbose, fprintf('* Time to filter random signals: %f sec.\n', t); end

if param.verbose, disp('Computing the SVD for eigenspace recovery...'); end
tic;
[basis, ~, ~] = svd(filt_sig, 'econ');
t = toc;
if param.verbose, fprintf('* Time to compute SVD: %f sec.\n', t); end

end


