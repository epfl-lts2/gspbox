function [ g ] = gsp_gabor_filterbank( G,k )
%GSP_GABOR_FILTERBANK returns a Gabor filterbank
%   Usage: g = gsp_gabor_filterbank( G,k );
%
%   Input parameters:
%       G       : Graph
%       k       : kernel
%
%   Output parameters:
%       g       : filterbank
%   
%   This function create a filterbank with the kernel k. Every filter is
%   centered in a different frequency
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_gabor_filterbank.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

% Author: Nathanael Perraudin
% Date  : 13 June 2014

if ~isfield(G,'E')
   error(['GSP_GABOR_FILTEBANK: You need first to compute the Fourier basis\n',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

Nf = length(G.E);

g = cell(Nf,1);

for ii=1:Nf
    g{ii} = @(x) k(x-G.E(ii));
end


end


