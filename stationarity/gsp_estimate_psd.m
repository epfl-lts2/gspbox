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
%    param.Nfilt  : Number of filters (default 50)
%    param.Nrandom : Number of random signal (default 10)
%    param.g0 : Initial window. Default:
%
%        g(x)  = exp(-( Nfilt^2 * x).^2/(Nfilt*lmax)^2) 
%
%   Additional parameters are availlable in the function
%   gsp_filter_analysis.
%   
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_estimate_psd.html

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


