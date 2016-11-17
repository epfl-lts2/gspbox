function [h,filtertype] = gsp_jtv_filter_array(G,g,filtertype)
%GSP_JTV_FILTER_ARRAY Convert ts/js filters to -array filters
%   Usage: [h,filtertype] = gsp_jtv_filter_array(G,g, filtertype)
%
%   Input parameters:
%         G          : Graph
%         g          : Cell array of time-vertex filters
%         filtertype : Filter domain (ts,js)
%   Output parameters:
%         h          : Cell array of graph filterbank
%         filtertype : Filter domain (ts-array,js-array)
%
%   Convert ts/js filters to -array filters
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_jtv_filter_array.php

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

% Author :  Francesco Grassi, Nathanael Perraudin
% Date : September 2016

if ~iscell(g)
   g = {g};
end

if ~gsp_check_filtertype(filtertype,{'ts','js'})
    error('Invalid filtertype.');
end

T  = G.jtv.T;
Nf = numel(g);

switch filtertype
    case 'ts'
        v = gsp_jtv_ta(G);
    case 'js'
        v = gsp_jtv_fa(G);
end

h  = cell(Nf,T);
for n=1:Nf
    for ii = 1:T
        h{n,ii} = @(x) g{n}(x,v(ii));
    end
end

filtertype = [filtertype '-array'];

end
