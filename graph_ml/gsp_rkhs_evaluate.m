function [ K ] = gsp_rkhs_evaluate( k,x,y, method )
%GSP_RKHS_EVALUATE Function that evaluate a kernel
%   Usage: [ K ] = gsp_rkhs_evaluate( k,x,y );
%  
%   Input parameters:
%       k       : kernel
%       x       : points (each feature vector is a column)
%       y       : points (each feature vector is a column) (default x)
%       method  : Method used to compute the kernel (default: 'mixed')
%
%   Output parameters:
%       K       : Evaluation of the kernel
%   
%   This function evaluate the kernel k between all point of x and y.
%   
%   Those different methods are availlable:
%   * 'copy'    : Copy both x and y an perform one shot evaluation of the
%                 kernel
%   * 'loop'    : Perform a double loop to evalute the kernl
%   * 'mixed'   : Use a single loop and a copy (this is the default)
%

% Author: Nathanael Perraudin
% Date  : 8 December 2014
    
if nargin < 4
    method = 'mixed';
end
    
if nargin <3
    y = x;
end

N = size(x,2);
M = size(y,2);

if iscell(k)
    K = zeros(N,M*length(k));
    for ii = 1:length(k)
        K(:,(ii-1)*M+(1:M)) = gsp_rkhs_evaluate( k{ii},x,y, method );
    end
    return
end


switch method
    case 'copy'
        K = k(repmat(x,M,1)',repmat(y,N,1));
    case 'loop'
        K = zeros(N,M);
        for ii = 1:N
            for jj = 1:M
               K(ii,jj) = k(x(:,ii),y(:,jj));
            end 
        end
    case 'mixed'
        K = zeros(N,M);
        for ii = 1:N
               K(ii,:) = k(repmat(x(:,ii),1,M),y);
        end
    otherwise
        error('GSP_RKHS_EVALUATE: Unknown method')
end




end

