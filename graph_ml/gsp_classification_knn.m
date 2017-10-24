function sol = gsp_classification_knn(G ,M, y )
%GSP_CLASSIFICATION_KNN Classification using knn
%   Usage: sol = gsp_classification_knn(G ,M, y );
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



%%

B = gsp_classification_matrix(y);

sol = gsp_regression_knn(G, M, B );

minf = min(y);
sol = gsp_matrix2label(sol,minf);



end