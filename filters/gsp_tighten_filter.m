function [ ftighten ] = gsp_tighten_filter(G, filters )
%GSP_TIGHTEN_FILTER Create a function that tighten a filterbank
%   Usage: ftighten = gsp_tighten_filter( filters );
%
%   Input parameters:
%       G           : Graph or maximum eigenvalue
%       filters     : Filters of the filterbank (cell array)
%   Ouput parameters:
%       ftighten    : Inline function
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using:
%
%       G = gsp_estimate_lmax(G);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_tighten_filter.php

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

% Author: Nathanael Perraudin, David Shuman
% Date  : 16 June 2014


if isstruct(G)
    if ~isfield(G,'lmax')
            warning('GSP_TIGHTEN_FILTER: has to compute lmax \n')
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end



[~,B] = gsp_filterbank_bounds([0,lmax],filters);

ftighten = @(x) f_tighten(filters,x,B);

end


function output = f_tighten(filters,x,B)
    output = B;

    for i=1:length(filters)
        output=output-(filters{i}(x)).^2;
    end
    output=real(sqrt(output));
end


