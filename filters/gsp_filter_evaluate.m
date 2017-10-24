function fd = gsp_filter_evaluate(filter, x)
%GSP_FILTER_EVALUATE Evaluate the filterbank
%   Usage: fd = gsp_filter_evaluate(filter, x)
%
%   Input parameters:
%       filter  : cell array of filter
%       x       : data
%   Output parameters:
%       fd      : response of the filters
%
%   This function applies all the filters in filter to the data x. Every
%   filter correspond to one column of the matrix fd.
%
%   See also: gsp_filter_analysis
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_filter_evaluate.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
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
% Date: 18 March 2014

if ~iscell(filter)
   filter = {filter};
end

Nf=numel(filter);
fd=zeros(length(x),Nf);
for k=1:Nf
    fd(:,k)=filter{k}(x);
end

end

