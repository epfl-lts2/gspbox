function [RMSE, N_obs] = gsp_rmse_mv(X, param)
%GSP_RMSE_MV Compute columnwise RMSE ignoring missing values (NaN's)
%   Usage: C = gsp_rmse_mv(X);
%          C = gsp_rmse_mv(X, param);
%          [C, N_obs] = gsp_rmse_mv(...);
%
%   Input parameters:
%       X           : Data matrix
%       param       : Parameter
%   Output parameters:
%       C           : RMSE
%       N_obs       : number of commonly observed values
%
%   The RMSE between two different columns will only take into account the
%   common observed values.
%
%   C(i,j) = xi-xj||/sqrt(N), where xi and xj only contain the elements
%   that are commonly observed in both and N is the number of these
%   elements. 
% 
%   If an element is equal to NaN, it is not taken into account.
%
%   N_obs gives the number of commonly observed values for all pairs of
%   columns.
%
%   param is an optional structure of argument containing the following
%   fields:
%    param.verbose*: Verbosity level of the function (default 0)
%
%   TODO: write fast implementation for case with no missing values!! like
%   corr of matlab does
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_rmse_mv.php

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


%
% Author: Vassilis Kalofolias, Nathanael Perraudin
% Date: Mar 2014
% Testing: test_rmse



if nargin < 2
    param = struct;
end

if isnumeric(param)
    scal = param;
    param = struct;
    param.verbose = scal;
end

%Parameters
if ~isfield(param,'verbose'), param.verbose = 0; end;

[N] = size(X, 2);

% keep the positions of the observed values
obs = not(isnan(X));


RMSE = zeros(N);

if param.verbose
    fprintf('RMSE_MV has begun:\n');
    print_cols = 10^(floor(log10(N)) - 1) * round(N / 10^floor(log10(N)));
    t = tic;
end

for i = 1:N-1
    % for all pairs i,j, j=1...N
    obs_ij = bsxfun(@and, obs(:,i), obs(:, i+1 : N));
    % number of common observed values for columns i and j (j = i...N)
    N_ij = full(sum(double(obs_ij)));
    % find all pairs of differences but keep only common elements before
    % summing
    
    X_ij = full(X(:, i+1 : N)); % this makes the next line faster;  
    X_ij = bsxfun(@plus, X_ij, -X(:, i));
    
    X_ij(not(obs_ij)) = 0;
    RMSE(i, i+1 : N) = sqrt(sum(X_ij.^2)./N_ij);
    
    if param.verbose
        if mod(i, print_cols) == 0
            fprintf('%d columns done after %d sec\n', i, toc(t));
        end
    end
end
    
% fill also the lower triangular part    
RMSE = RMSE + RMSE';

if param.verbose

	fprintf('RMSE_MV has finished after %d sec:\n', toc(t));

end

if nargout>1
    N_obs = double(obs)'*double(obs);
end

end



