function [sol, infos] = gsp_wiener_inpainting(G,y, M, psd, psd_noise, param)
%GSP_WIENER_INPAINTING Solve wiener in-painting problem
%   Usage:  sol = gsp_wiener_inpainting(G,y, M, psd, psd_noise)
%           sol = gsp_wiener_inpainting(G,y, M, psd, psd_noise, param)
%           [sol, infos] = gsp_wiener_inpainting(...)
%
%   Input parameters:
%         G          : Graph (GSP structure)
%         y          : Measurements (column vector)
%         M          : Mask (vector)
%         psd        : PSD filter (anonymous function)
%         psd_noise  : PSD filter of the noise or single number
%         param      : Optional optimization parameters
%   Output parameters:
%         sol        : Solution
%         infos      : Convergence informations
%
%   This function solves the following wiener optimization problem:
%
%        argmin_x || M x - y ||_2^2 + || w(L) x ||_2^2 
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
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_wiener_inpainting.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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

Mop =@(x) bsxfun(@times,M,x);

if isnumeric(psd_noise) && sum(abs(psd_noise))==0
    ffid.eval = @(x) eps;
    ffid.prox = @(x,T) x-Mop(x)+Mop(y);
    wl = @(x) 1./(psd(x)+eps);
    fprox = @(T) @(x) psd(x)./(2*T + psd(x) + eps);  
else
    ffid.grad = @(x) 2*Mop(Mop(x)-y);
    ffid.eval = @(x) norm(Mop(x)-y,'fro')^2;
    ffid.beta = 2;
    
    if isnumeric(psd_noise)
        wl = @(x) psd_noise./(psd(x)+eps);
        fprox = @(T) @(x) psd(x)./(psd(x)+2*T*psd_noise + eps);  
    else
        wl = @(x) psd_noise(x)./(psd(x)+eps);
        fprox = @(T) @(x) psd(x)./(psd(x)+2*T*psd_noise(x) + eps);
    end
        
end

    
fprior.prox = @(x,T) gsp_filter_analysis(G,fprox(T),x, param);
fprior.eval = @(x) 0.5*norm(gsp_filter_analysis(G,wl,x,param),'fro')^2;


[sol, infos] = solvep(y,{ffid,fprior},param);
% [sol, infos] = gsp_wiener_optimization(G, y, ffid, psd, psd_noise, param);

end
