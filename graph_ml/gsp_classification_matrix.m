function B = gsp_classification_matrix(f)
%GSP_CLASSIFICATION_MATRIX Create the classification matrix from label f
%   Usage: B = gsp_classification_matrix(f);
%
%   Input parameters:
%       f   : Labels
%   Output parameters:
%       B   : Classification matrix
%
%   See also: gsp_matrix2label gsp_classification_tik gsp_classification_tv
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graph_ml/gsp_classification_matrix.html

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
% Date  : 24 July 2015

minf = min(f);
maxf = max(f);
d = minf:maxf;
B = zeros(numel(f),numel(d));
for ii = 1:numel(d)
    B(:,ii) = double(f==d(ii));
end

B = double(B);

end
