function psd = gsp_estimate_psd(G, x, param)
%GSP_ESTIMATE_PSD Estimation of the Power Spectrum Density
%   Usage:  psd = gsp_estimate_psd(G, x)
%           psd = gsp_estimate_psd(G, x, param)
%
%   Input parameters:
%         G          : Graph
%         x          : Signal(s) (vector column)
%         param      : Optional parameters
%   Output parameters:
%         psd        : PSD filter
%
%   This function estimate the PSD from a signal using the method presented
%   in "Stationary signal processing on graphs"
%
%   Additional parameters
%   ---------------------
%  
%   * *param.Nfilt*  : Number of filters (default 50)
%   * *param.Nrandom* : Number of random signal (default 10)
%   * *param.g0* : Initial window. Default:
%
%     .. g(x)  = exp(-( Nfilt^2 * x).^2/(Nfilt*lmax)^2) 
%
%     .. math:: g(x)  = e^{-\frac{N_f^4 x^2}{(N_f + 1)^2\lambda_{\max}^2}}
%
%   Additional parameters are availlable in the function
%   gsp_filter_analysis.
%   
%   References: perraudin2016stationary

% Author : Nathanael Perraudin
% Date: 6 January 2016


%% Handling optional parameters
if nargin<3
    param = struct;
end

if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);

    warning(['GSP_PSD_ESTIMATION: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
end

if ~isfield(param,'Nfilt'), param.Nfilt = 50; end
if ~isfield(param,'Nrandom'), param.Nrandom = max(10,size(x,2)); end
if ~isfield(param,'g0'), 
    sigma = sqrt(2*G.lmax/param.Nfilt^2 * (param.Nfilt + 1));
    param.g0 = @(x) exp(-x.^2/sigma^2); 

end
%% Design the frame

[ g , mu ] = gsp_design_translates(G, param.g0, param.Nfilt);
%[ g , mu ] = gsp_design_itersine(G,Nfilt );

%% Perform the filtering

x_filt = gsp_vec2mat( gsp_filter(G,g,x,param), param.Nfilt);
%% estimate the points
n2_x_filt = sum(abs(x_filt).^2,1);
mu_y = reshape(mean(n2_x_filt,3),[],1);

%% Estimate the energy of the window

if gsp_check_fourier(G)
   mu_y2 = sum(abs(gsp_filter_evaluate(g,G.e)).^2,1)';
else
    w = randn(size(x,1),param.Nrandom);
    x_filt2 = gsp_vec2mat( gsp_filter(G,g,w,param), param.Nfilt);
    n2_x_filt2 = sum(abs(x_filt2).^2,1);
    mu_y2 = reshape(mean(n2_x_filt2,3),[],1);
end


%% Interpolate to obtain a nice filter.

psd = @(s) max(spline(mu,mu_y./mu_y2,s),0);

end

