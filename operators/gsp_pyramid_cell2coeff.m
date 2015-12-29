function coeff = gsp_pyramid_cell2coeff(ca,pe)
%GSP_PYRAMID_CELL2COEFF Cell array to vector transform for the pyramid
%   Usage: coeff = gsp_pyramid_cell2coeff(ca,pe);
%
%   Input parameters:
%       ca      : Cell array with the coarse approximation at each level
%       pe      : Cell array with the prediction errors at each level
%   Output parameters:
%       coeff   : Vector of coefficient
%
%   This function compress the cell array ca and pe into a single
%   vector of coefficients. It keeps  the smaller coarse approximation and
%   the prediction errors. 
%
%   Example:
%
%           [ca,pe] = gsp_pyramid_analysis(Gs, f);
%           coeff = gsp_pyramid_cell2coeff(ca,pe);
%
%   See also: gsp_pyramid_analysis gsp_pyramid_synthesis
%             gsp_graph_multiresolution
%
%   Demo: gsp_demo_pyramid
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_pyramid_cell2coeff.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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
% Date  : 5 Aout 2014
% Testing : test_pyramid

Nl = length(ca) - 1;

N = 0;
for ii = 1 : (Nl+1)
    N = N+length(ca{ii});
end

Nv = size(ca{Nl}, 2);
coeff = zeros(N,Nv);

Nt = length(ca{Nl+1});
coeff(1:Nt, :) = ca{Nl+1};

ind = Nt+1;
for ii = 1 : Nl
    Nt = length(ca{Nl+1-ii});
    coeff(ind:(ind+Nt-1), :) = pe{Nl+1-ii};
    ind = ind + Nt;
end

if (ind-1) ~= N 
    error('Something is wrong here: contact the gspbox team')
end

end


