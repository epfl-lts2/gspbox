function [ fa ] = gsp_cfa( T,fs )
%GSP_CFA Create frequency axis
%   Usage: fa = gsp_cfa(N);
%          fa = gsp_cfa(N,fs);
%
%   Input parameters:
%       N   : Number of samples
%       fs  : Sampling frequency (default 1)
%   Ouput parameters:
%       fa  : Frequency axis
%   
%   This function create a normalized frequency axis corresponding to the
%   frequency of the function fft.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_cfa.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

if nargin<2
    fs = 1;
end

if mod(T,2)
    fa = linspace(-0.5+1/(2*T),0.5-1/(2*T),T);
else
    fa = linspace(-0.5,0.5-1/T,T);
end

fa = fa*fs;

fa = ifftshift(fa);
fa = fa(:);

end

