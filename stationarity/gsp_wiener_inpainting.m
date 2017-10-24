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
%     .. argmin_x || M x - y ||_2^2 + || w(L) x ||_2^2 
%
%     .. math:: arg\min_x \| M x - y \|_2^2 + \| w(L) x \|_2^2 
%
%   Please refer to the reference for more information about this problem.
%   This function requires the UNLocBox to work.
%
%   Please refer to the function gsp_filter_analysis and solvep to know how
%   *param* can be set.
%
%   References: perraudin2016stationary

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