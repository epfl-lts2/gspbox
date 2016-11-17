function [Wsym,Wasym] = gsp_decompose_asymmatrix(W)
%GSP_DECOMPOSE_ASYMMATRIX Decompose a matrix in symmetric and asymmetric part
%   Usage: [Wsym,Wasym] = gsp_decompose_asymmatrix(W)
%
%   Input parameters:
%       W        : Asymmetric matrix
%   Output parameters:
%       Wsym    : Symmetric part of W
%       Wasym   : Asymmetric part of W
%
%   Decompose a matrix in symmetric and asymmetric part
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_decompose_asymmatrix.php

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

% Author: Francesco Grassi
% Date  : July 2016

Wasym = W-W.';

Wasym(Wasym<0)=0;

Wsym = W-Wasym;

Wsym = (Wsym+Wsym')/2;

end
