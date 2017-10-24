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
%   |gsp_knn_classify_graph| .
%
%   See also: gsp_classification_tik gsp_classification_tv

% Author: Nathanael Perraudin
% Date  : 24 July 2015
% Testing: test_graph_ml

    solt = repmat(G.de.^(-1),1,size(y,2)) .* (G.We * y(logical(M),:) );
    sol = y;
    sol(logical(1-M),:) = solt;
end