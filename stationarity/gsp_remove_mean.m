function [ X, mX ] = gsp_remove_mean( X, dim )
%GSP_REMOVE_MEAN Remove the mean of the signal X
%   Usage:  X  = gsp_remove_mean( X );
%           X  = gsp_remove_mean( X, dim );
%           Ç X, mX ]  = gsp_remove_mean( ... );
%
%   Input parameters:
%       X       : data
%       dim     : dimension (default 1)
%   Output parameters:
%       X       : centered data
%       mX      : mean
%   
%   This function remove efficiently the mean of X across line or column of
%   the data.
%

% Author: Nathanael Perraudin
% Date  : 24 July 2016
% Testing: test_gsp_remove_mean



if nargin<2
    dim = 1;
end

mX = mean(X,dim);

X = bsxfun(@minus,X,mX); 


end

