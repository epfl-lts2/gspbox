function gsp_plot_jft(G,Xhat,param)
%GSP_PLOT_JFT  Plot the magnitude squared of the joint Fourier transform matrix Xhat
%   Usage:  gsp_plot_jft(G,X);
%           gsp_plot_jft(G,X,param);
%
%   Input parameters:
%       G       : Time-Vertex graph structure
%       Xhat    : Joint Time-Vertex Fourier Coefficients Matrix
%       param   : Structure of optional parameters
%   Output parameters:
%       none
%
%   Additional parameters
%   ---------------------
%
%    param.dim       : '2d' for imagesc and '3d' for surf (default 3d)
%    param.logscale  : Use log-scale to visualize the JFT (default 0)
%    param.dB        : Range of plotting in dB for the logscale
%    param.fftshift  : Put zero frequency at center (default 1)
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/plotting/gsp_plot_jft.html

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

% Author: Francesco Grassi, Nathanael

if nargin<3
    param=struct;
end

if ~isfield(param,'dim');      param.dim='3d';end
if ~isfield(param,'logscale'); param.logscale=0; end
if ~isfield(param,'dB');       param.dB = Inf; end
if ~isfield(param,'fftshift'); param.fftshift = 1; end


if isfield(G,'e')
    lambda = G.e;
else
    if isfield(G,'lmax')
        lambda = linspace(0,G.lmax,G.N);
    else
        lambda = linspace(0,G.N,G.N);
    end
end




if param.logscale
    Xhat = 20*log(abs(Xhat)+1);
    maxXhat = max(Xhat(:));
    Xhat(Xhat < (maxXhat-param.dB)) = maxXhat-param.dB;
else
    Xhat=abs(Xhat.^2);
end

if param.fftshift
    omega = gsp_jtv_fa( G,1 );
    Xhat  = fftshift(Xhat,2);
else
    omega = gsp_jtv_fa( G,0 );
end


switch param.dim
    case {2,'2','2d'}
        imagesc(omega,lambda,Xhat);
        set(gca,'YTickLabel','')
        axis xy
        colorbar
        
    case {3,'3','3d'}
        surf(omega,lambda,Xhat,'linestyle','none');
        view([0 90])
        xlim([min(omega) max(omega)])
        ylim([min(lambda) max(lambda)])
        
    otherwise
        error('unkown plot method');
end


