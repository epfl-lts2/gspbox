function [ ft ] = gsp_translate(G, f, i)
%GSP_TRANSLATE Generalized translation of the signal f to the node i
%   Usage: ft = gsp_translate(G, f, i);
%
%   Input parameters
%       G   : Graph
%       f   : Signal (column)
%       i   : Indices of vertex (int)
%   Output parameters
%       ft  : translate signal
%
%   This function translate the column vector f onto the node i. If f*
%   is a matrix, the translation will be done to each column.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_translate.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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


fhat=gsp_gft(G,f);
nt = size(f,2);

ft = sqrt(G.N)*gsp_igft(G,fhat .* ...
    repmat(transpose(G.U(i,:)),1,nt));


end
