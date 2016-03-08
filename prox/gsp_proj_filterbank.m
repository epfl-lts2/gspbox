function [sol, infos] = gsp_proj_filterbank(x, ~ , G, W, y, param)
%GSP_PROJ_FILTERBANK Projection onto the synthesis coefficients
%   Usage:  sol = gsp_proj_filterbank(x, 0, G, W, y, param);
%           sol = gsp_proj_filterbank(x, 0, G, W, y);
%           [sol, info] = gsp_proj_filterbank(...)
%
%   Input parameters:
%         x     : Input signal
%         G     : Graph structure
%         W     : Filterbank (cell array of functions)
%         y     : Measurements
%         param : Structure of optional parameters
%   Output parameters
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   GSP_PROJ_FILTERBANK(x, gamma, G, W, param) can solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2  such that W^* x =y 
%
%   Where W is the linear analysis operator associated with the
%   filterbank. 
%
%   param is a Matlab structure containing the following fields:
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%    param.weights : weights for a weighted L2-norm (default = 1)
%
%   info is a Matlab structure containing the following fields:
%
%    infos.algo : Algorithm used
%    infos.iter : Number of iteration
%    infos.time : Time of exectution of the function in sec.
%    infos.final_eval : Final evaluation of the function
%    infos.crit : Stopping critterion used 
%
%
%   See also:  gsp_solve_l1 gsp_proj_b2_filterbank 
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/prox/gsp_proj_filterbank.php

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


% Author: Nathanael Perraudin
% Date: 25 March 2014
% Testing: test_gsp_proj_filerbank




if nargin < 5
    error('GSP_PROJ_FILTERBANK: Not enought input arguments');
end

if nargin < 6, param=struct; end

if ~isfield(param, 'eps'), param.eps = 1e-8; end


if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    warning(['GSP_PROJ_FILTERBANK: To be more efficient you should run: ',...
        'G = gsp_estimate_lmax(G); before using this proximal operator.']);
end

t = tic;

if size(y,2)==1 && size(x,2)>1
    y = repmat(y,1,size(x,2));
end
Wd = gsp_design_can_dual(W,param.eps);
sol =  x - gsp_filter_analysis(G,Wd,(gsp_filter_synthesis(G,W,x)-y));


infos.iter = 1;
infos.time = toc(t);
infos.algo = mfilename;
infos.final_eval = 0;
infos.crit = 'exact';


end



