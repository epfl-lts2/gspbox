function sol = gsp_jtv_wiener_inpainting(G, y, M, psd, psdnoise, param)
%GSP_JTV_WIENER_INPAINTING Solve wiener in-painting problem for vertex time process
%   Usage:  sol = gsp_jtv_wiener_inpainting(G, y, M, psd, psdnoise)
%           sol = gsp_jtv_wiener_inpainting(G, y, M, psd, psdnoise, param)
%           [sol, infos] = gsp_jtv_wiener_inpainting(...)
%
%   Input parameters:
%         G          : Graph (GSP structure)
%         y          : Measurements
%         M          : Mask
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
% Testing : test_gsp_wiener_jft_predition


if nargin<6
    param = struct;
end

if ~isfield(param, 'accelerate')
    param.accelerate = 1;
end

if param.accelerate && ~isnumeric(psd) && isnumeric(psdnoise) && gsp_check_fourier(G)
    if numel(psd)>1
        psd = gsp_filter_evaluate(psd,G.e);
    else
        psd = gsp_jtv_filter_evaluate(psd,'js',G.e,gsp_jtv_fa(G));
    end
    
end

if isnumeric(psd) && isnumeric(psdnoise)
    if sum(abs(psdnoise(:)))>0
        wiener_filt =@(T) psd./(2*psdnoise*T+psd+eps);
        ffid.eval = @(x) norm(M.*x-y,'fro')^2;
        ffid.grad = @(x) 2*M.*(M.*x-y);
        ffid.beta = 2;
        weights = psdnoise./(psd+eps);
        
    else
        wiener_filt =@(T) psd./(2*T+psd+eps);
        ffid.eval = @(x) eps;
        ffid.prox = @(x,T) x-M.*x+M.*y;
        weights = 1./(psd+eps);
        
    end
    
    fprior.prox = @(x,T) gsp_ijft(G,wiener_filt(T).*gsp_jft(G,x));
    fprior.eval = @(x) norm(weights.*gsp_jft(G,x),'fro')^2;
else
    
    if numel(psd)>1
        ft = 'js-array';
        if isnumeric(psdnoise) && sum(abs(psdnoise))==0
            ffid.eval = @(x) eps;
            ffid.prox = @(x,T) x-M.*x+M.*y;
            
            %cell2mat(cellfun(@(h) h(x),h(ii,:),'uniformoutput',0));
            %wl = @(lambda,omega) 1./(psd(lambda,omega)+eps);
            wl = apply2array(psd, @(x) 1./(x+eps));
            %fprox = @(T) @(lambda,omega) psd(lambda,omega)./(2*T + psd(lambda,omega) + eps);
            fprox = @(T) apply2array(psd, @(x) x./(2*T + x + eps));
        else
            ffid.grad = @(x) 2*M.*(M.*x-y);
            ffid.eval = @(x) norm(M.*x-y,'fro')^2;
            ffid.beta = 2;
            
            if isnumeric(psdnoise)
                % wl = @(lambda,omega) psdnoise./(psd(lambda,omega)+eps);
                wl = apply2array(psd, @(x) psdnoise./(x+eps));
                % fprox = @(T) @(lambda,omega) psd(lambda,omega)./(psd(lambda,omega)+2*T*psdnoise + eps);
                fprox = @(T) apply2array(psd,@(x)  x./(2*T*psdnoise + x + eps));
            else
                error('Not implemented yet!')
            end
            
        end
        
    else
        ft = 'js';
        
        if isnumeric(psdnoise) && sum(abs(psdnoise))==0
            ffid.eval = @(x) eps;
            ffid.prox = @(x,T) x-M.*x+M.*y;
            wl = @(lambda,omega) 1./(psd(lambda,omega)+eps);
            fprox = @(T) @(lambda,omega) psd(lambda,omega)./(2*T + psd(lambda,omega) + eps);
        else
            ffid.grad = @(x) 2*M.*(M.*x-y);
            ffid.eval = @(x) norm(M.*x-y,'fro')^2;
            ffid.beta = 2;
            
            if isnumeric(psdnoise)
                wl = @(lambda,omega) psdnoise./(psd(lambda,omega)+eps);
                fprox = @(T) @(lambda,omega) psd(lambda,omega)./(psd(lambda,omega)+2*T*psdnoise + eps);
            else
                wl = @(lambda,omega) psdnoise(lambda,omega)./(psd(lambda,omega)+eps);
                fprox = @(T) @(lambda,omega) psd(lambda,omega)./(psd(lambda,omega)+2*T*psdnoise(lambda,omega) + eps);
            end
            
        end
    end
    fprior.prox = @(x,T) gsp_jtv_filter_analysis(G,fprox(T),ft,x, param);
    fprior.eval = @(x) 0.5*norm(gsp_jtv_filter_analysis(G,wl,ft,x,param),'fro')^2;
end

sol = solvep(y,{ffid,fprior},param);


end


function g = apply2array(garray,fun)

g = cell(size(garray));
for ii = 1:numel(garray)
    g{ii} = @(x) fun(garray{ii}(x));
end

end