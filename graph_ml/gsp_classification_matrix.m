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