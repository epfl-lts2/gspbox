function [Wsym,Wasym] = gsp_decompose_asymmatrix(W)
%GSP_DECOMPOSE_ASYMMATRIX Decompose a matrix in symmetric and asymmetric part
%   Usage: [Wsym,Wasym] = gsp_decompose_asymmatrix(W)
%
%   Input parameters:
%       W        : Asymmetric matrix
%   Output parameters:
%       Wsym    : Symmetric part of W
%       Wasym   : Asymmetric part of W
%
%   Decompose a matrix in symmetric and asymmetric part
%

% Author: Francesco Grassi
% Date  : July 2016

Wasym = W-W.';

Wasym(Wasym<0)=0;

Wsym = W-Wasym;

Wsym = (Wsym+Wsym')/2;

end