function C = gsp_stationarity_cov(X)
%GSP_STATIONARITY_COV Covariance matrix from graph stationary data
%   Usage:  C = gsp_stationarity_cov(X)
%
%   Input parameters:
%         X          : Data (M x N matrix) 
%   Output parameters:
%         C          : Covariance matrix (M x M)
%
%   This function estimates the covariance from the data. Every sample has
%   the same expected average.
%
%   References: perraudin2016stationary

% Author : Nathanael Perraudin
% Date: 6 January 2016

Nso = size(X,2);
X = X - mean(X(:));
C = (X*X')/Nso;


end