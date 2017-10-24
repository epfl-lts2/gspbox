function [ s ] = gsp_evaluate_can_dual( g,val,tol )
%GSP_EVALUATE_CAN_DUAL Evaluate the canonical dual filterbank
%   Usage: hcoeff = gsp_evaluate_can_dual( g,val )
%
%   Inputs parameters:
%       g       : cell array of filters
%       val     : column vectors of values
%       tol     : tolerance
%
%   Ouputs parameters:
%       s       : Matrix of value
%
%   This function computes the value of the canonical dual of a filterbank
%   *g* at the point specified in *val*. The function returns a matrix. Each
%   column is the output of one dual filter.
%
%   See also: gsp_design_can_dual

% Author: Nathanael Perraudin
% Date  : 14 June 2014
% Testing: test_dual

if nargin<3
    tol = 1e-10;
end


% TODO: size should be improved
N = length(val);

% Compute coefficient of g
gcoeff = gsp_filter_evaluate(g,val)';

% Compute coefficient of h
% s = zeros(N,M);
% for ii = 1:N
%     s(ii,:) =  pinv(gcoeff(ii,:)'); 
% end

s = arrayfun(@(x) pinv(gcoeff(:, x),tol), 1 : N, 'UniformOutput', false);
s = cell2mat(s');
end

