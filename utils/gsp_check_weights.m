function [ a ] = gsp_check_weights( W )
%GSP_CHECK_WEIGHTS Check a weight matrix
%   Usage: a = gsp_check_weights( W );
%
%   Input parameters:
%       W       : Weight matrix
%   Output parameters:
%       a       : Warning code
%   
%   This function performs various test on the weight matrix $W$. It
%   returns:
%
%       * 0     : Everything is ok
%       * 1     : The martrix contains inf values
%       * 2     : The diagonal is not  0
%       * 3     : The matrix is not square
%       * 4     : The matrix contains nan values
%

% Author: Nathanael Perraudin
% Date  : 12 june 2014

a = 0;


if sum(sum(isinf(W)))
    warning(['GSP_CHECK_WEIGHTS: There is infinite value',...
             ' in the weight matrix']);
    a = 1;
end

if sum(abs(diag(W)))
    disp(['GSP_CHECK_WEIGHTS: The diagonal',...
             ' of the weight matrix is not 0!']);
    a = 2;
end

if size(W,1) ~= size(W,2)
    warning('GSP_CHECK_WEIGHTS: The weight matrix is not square!');
    a = 3;
end

if sum(sum(isnan(W)))
    warning(['GSP_CHECK_WEIGHTS: There is infinite value',...
             ' in the weight matrix']);
    a = 4;
end

end

