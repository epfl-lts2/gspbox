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
%     .. argmin_x f(x) + || w(L) x ||_2^2 
%
%     .. math:: arg\min_x f(x) + \| w(L) x \|_2^2 
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