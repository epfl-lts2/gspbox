function [ C ] = gsp_gwft2(G,f,k, param )
%GSP_GWFT Graph windowed Fourier transform
%   Usage:  G = gsp_gwft(G,f,k, param );
%           G = gsp_gwft(G,f,k);
%
%   Input parameters:
%         G     : Graph
%         f     : Graph signal
%         k     : kernel
%         param : Structure of optional parameter
%   Output parameters:
%         C     : Coefficient.
%
%   This function compute the graph windowed Fourier transform of a signal
%   f with the window g. The function returns a matrix of size N^2*N.
%   This function take k as a kernel and translate it a all frequencies.
%
%   param a Matlab structure containing the following fields:
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%     By default, it is 1.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_gwft2.php

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

% Author : Nathanael Perraudin
% Testing: test_gwft


% Optional input arguments
if nargin<4, param=struct; end

if ~isfield(param, 'verbose'), param.verbose=1 ; end


if ~isfield(G,'E')
   error(['GSP_GABOR_GWFT2: You need first to compute the Fourier basis\n',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

g = gsp_gabor_filterbank( G,k );


C = gsp_filter_analysis(G,g,f);

C = transpose(gsp_vec2mat(C,G.N));
    
end


