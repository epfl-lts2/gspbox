function sol = gsp_regression_knn(G,M,y)
%GSP_REGRESSION_KNN Regression using knn
%   Usage: sol = gsp_regression_knn(G ,M, y );
%
%   Input parameters:
%       G   : Graph
%       M   : Mask (to determine with label is known)
%       y   : label (total size of the problem)
%
%   Output parameters:
%       sol : Solution of the problem
%
%   Warning: In order to use this function, you have to use a special
%   graph. This graph can be computed with the function:
%   GSP_KNN_CLASSIFY_GRAPH .
%
%   See also: gsp_classification_tik gsp_classification_tv
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graph_ml/gsp_regression_knn.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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
% Testing: test_graph_ml

    solt = repmat(G.de.^(-1),1,size(y,2)) .* (G.We * y(logical(M),:) );
    sol = y;
    sol(logical(1-M),:) = solt;
end
