function W = gsp_symmetrize(W, type)
%GSP_SYMMETRIZE symmetrize a matrix
%   Usage:  W = gsp_symmetrize(W)
%           W = gsp_symmetrize(W, type)
%
%   Input parameters:
%       W       : square matrix
%       type    : type of symmetrization (default 'full')
%   Output parameters:
%       W       : symmetrized matrix
%
%   The available symmetrization types are:
%    'average' : average of W and W^T (default)
%    'full'    : copy the missing entries
%    'none'    : nothing is done (the matrix might stay unsymmetric!)
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_symmetrize.php

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
% Date  : 17 Janvier 2015

if nargin < 2
    type = 'full';
end

switch type
    case 'average'
        W = (W+W.')/2;
    case 'full'
        A = W>0;
        M = logical(A - (A' & A));

        W = W + M'.*W';
 %       Wt = W';
 %       W(M') = Wt(M');
        
    case 'none'
        return
    otherwise
        error('GSP_SYMMETRIZE: Unknown type')
end

end

