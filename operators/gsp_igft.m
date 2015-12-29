function [f] = gsp_igft(G,f_hat)
%GSP_IGFT  Inverse graph Fourier transform
%   Usage:  f = gsp_igft(G,f_hat);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         f_hat      : Signal
%   Output parameters:
%         f          : Inverse graph Fourier transform of f_hat*
%
%   'gsp_igft(G,f_hat)' computes a graph Fourier transform of the signal
%   f_hat with respect to the Fourier basis of the graph G: G.U.
%   Alternatively, one can provide directly the Fourier basis instead of
%   the graph G. 
%
%      f = U * f_hat 
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
%       f_hat = zeros(N,1);
%       f_hat(5) = 1;
%       f = gsp_igft(G,f_hat);
%       gsp_plot_signal(G,f);  
%
%   See also: gsp_gft, gsp_gwft, gsp_compute_fourier_basis
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_igft.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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
       error(['GSP_IGFT: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    U = G.U;
else
    U = G;
end

f=U*f_hat;

