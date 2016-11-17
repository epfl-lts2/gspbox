function G = gsp_incidence(G,inc_type)
%GSP_INCIDENCE Compute an incidence matrix
%
%   Usage: G = gsp_incidence(G)
%              gsp_incidence(G,inc_type)
%
%   Input parameters:
%       G        : Graph structure
%       inc_type : Type of incidence matrix: 'weighted' or 'binary'
%   Output parameters:
%       G   : Graph structure
%
%   'gsp_incidence' compute the incidence matrix B from the gradient and add it to graph structure
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_incidence.php

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

if nargin < 2
    inc_type = 'weighted';
end


if ~isfield(G,'Diff'); G = gsp_adj2vec(G); end;

switch inc_type
    case 'weighted'
        B = abs(G.Diff);
    case 'binary'
        B = abs(G.Diff)>0;
    otherwise
        error('Unknown incidence matrix.')
end

G.B = B;

end
