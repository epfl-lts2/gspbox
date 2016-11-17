function bool = gsp_check_jtv(G)
%GSP_CHECK_JTV Check if G is a JTV graph
%   Usage:  bool = gsp_check_jtv(G):
%
%   Input parameters:
%       G           : Graph structure
%   Output parameters:
%       bool        : boolean
%
%   This function check if the structure G is a valid Joint Time-Vertex
%   Graph structure
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_check_jtv.php

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

bool = 0;

if ~isstruct(G)
    error('Input is not a valid Graph structure')
end

if ~isfield(G.jtv,'T')
    return
end

if or(~isfield(G.jtv,'T'),~isfield(G.jtv,'fs'))
    return
end

bool = 1;
