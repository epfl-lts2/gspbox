function [sol, infos] = gsp_wiener_optimization(G, x0, f, psd, psd_noise, param)
%GSP_WIENER_OPTIMIZATION Solve wiener optimization problem
%   Usage:  sol = gsp_wiener_optimization(G, x0, ffid, psd, psd_noise)
%           sol = gsp_wiener_optimization(G, x0, ffid, psd, psd_noise, param)
%           [sol, infos] = gsp_wiener_optimization(...)
%
%   Input parameters:
%         G          : Graph (GSP structure)
%         x0         : Starting point (column vector)
%         f          : Fidelity term - UNLocBox structure
%         psd        : PSD filter (anonymous function)
%         psd_noise  : PSD filter of the noise or single number
%         param      : Optional optimization parameters
%   Output parameters:
%         sol        : Solution
%         infos      : Convergence informations
%
%   This function solves the following wiener optimization problem:
%
%        argmin_x f(x) + || w(L) x ||_2^2 
%
%   Please refer to the reference for more information about this problem.
%   This function requires the UNLocBox to work.
%
%   Please refer to the function gsp_filter_analysis and solvep to know how
%   param can be set.
%   
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/stationarity/gsp_wiener_optimization.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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

% Author : Nathanael Perraudin
% Date: 6 January 2016


if nargin<6
    param = struct;
end

if isnumeric(psd_noise)
    if sum(abs(psd_noise(:)))==0
        error('This function cannot solve this problem')
    end
    if sum(abs(psd_noise(:)))<1e-10
        warning('This function can prabaly not solve this case');
    end
    wl = @(x) psd_noise./(psd(x)+eps);
    fprox = @(T) @(x) psd(x)./(psd(x)+2*T*psd_noise + eps);           
else
    wl = @(x) psd_noise(x)./(psd(x)+eps);
    fprox = @(T) @(x) psd(x)./(psd(x)+2*T*psd_noise(x) + eps);

end

%fprox = @(x) 1./(wl(x)+1);

% In order to be faster
param.stopping_criterion = 'rel_norm_obj';


% Wiener term 
fwiener.prox = @(x,T) gsp_filter_analysis(G,fprox(T),x, param);
fwiener.eval = @(x) 0.5*norm(gsp_filter_analysis(G,wl,x,param),'fro')^2;

% Call the solver
[sol , infos ] = solvep(x0,{f,fwiener},param);

end
