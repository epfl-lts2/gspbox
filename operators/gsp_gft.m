function [f_hat]=gsp_gft(G,f)
%GSP_GFT  Graph Fourier transform
%   Usage:  f_hat=gsp_gft(G,f);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         f          : f (signal)
%   Output parameters:
%         f_hat      : Graph Fourier transform of f*
%
%   'gsp_gft(G,f)' computes a graph Fourier transform of the signal f
%   with respect to the Fourier basis of the graph G: G.U. Alternatively,
%   one can provide directly the Fourier basis instead of the graph G.
%
%      f_hat = U' * f 
%
%   To compute the Fourier basis of a graph G, you can use the function:
%
%           G = gsp_compute_fourier_basis(G);
%
%   Example:
%
%       N = 30;
%       G = gsp_sensor(N);
%       G = gsp_compute_fourier_basis(G);
%       f = sin((1:N)'*2*pi/N);
%       fhat = gsp_gft(G,f);
%       gsp_plot_signal_spectral(G,fhat);  
%
%   See also: gsp_igft, gsp_gwft, gsp_compute_fourier_basis
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_gft.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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

% Author : Nathanael Perraudin, David I Shuma
% Date: 


if ~isnumeric(G);
    if ~isfield(G,'U')
       error(['GSP_GFT: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    U = G.U;
else
    U = G;
end

f_hat=U'*f;

