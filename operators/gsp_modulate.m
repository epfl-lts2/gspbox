function [ fm ] = gsp_modulate( G,f,k )
%GSP_MODULATE Tranlate the signal f to the node i
%   Usage: fm = gsp_modulate( G,f,k );
%
%   Input parameters
%       G   : Graph
%       f   : Signal (column)
%       k   : Indices of frequencies (int)
%   Output parameters
%       fm  : Modulated signal
%
%   This function modulate the column vector f onto the node i. If f is a
%   matrix, the modulation will be applicated to each column.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_modulate.php

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

% Author: Nathanael Perraudin
% Date  : 09.12.2013

nt = size(f,2);

fm = sqrt(G.N)*repmat(f,1,nt).*repmat(G.U(:,k+1),1,nt);


end

