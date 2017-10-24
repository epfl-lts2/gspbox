function [ f ] = gsp_matrix2label( B, minf )
%GSP_MATRIY2LABEL Reconstruct labels from matrix
%   Usage:  f = gsp_matrix2label( B );
%           f = gsp_matrix2label( B, minf );
%   
%   Input parameters:
%       B   : Classification matrix
%       minf: smallest integer (default 0)
%   Output parameters:
%       f   : Labels
%
%   See also: gsp_classification_matrix gsp_classification_tik gsp_classification_tv

% Author: Nathanael Perraudin
% Date  : 24 July 2015

if nargin<2
    minf = 0;
end

[~,f] = max(B,[],2);
f = f-1+minf;

end

