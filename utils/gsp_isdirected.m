function bool = gsp_isdirected(G)
%GSP_ISDIRECTED Check is the graph is directed
%   Usage: bool = gsp_isdirected(G);
%          bool = gsp_isdirected(W);
%
%   Input parameters
%       G       : Graph structure or square matrix
%   Output parameters
%       bool    : Boolean
%
%   This function test if the graph is directed. Alternatively, you can
%   give a square matrix and it tests if it is symetric. The function
%   returns 0 if the matrix is symetric and 1 otherwise!
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_isdirected.php

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
% Date  : 13 August 2014

if isstruct(G)
    W = G.W;
else
    W = G;
end

bool = sum(sum( abs(W - transpose(W))> eps(10) ))>0;

end

