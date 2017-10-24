function [x,y] = gsp_point2dcdf(v)
%GSP_POINT2DCDF points to discrete continuous density function
%   Usage : [x,y] = gsp_point2dcdf(v);
%
%   Input parameters:
%       v   : vector of sample (1 dimention)
%   Output parameters:
%       x   : coordinate along x
%       y   : coordinate along y
%
%   This function compute a discrete continuous density function from a
%   sample list. 
%   
%   Example:::
%
%       gsp_point2dcdf([0 3 4 2 2 0 1 0 1 1 2])
%

N = numel(v);
tol = 1e-10;
v = sort(v);
[x, inds] = unique(round(v(:)*1/tol)*tol);

y = (inds-1)/(N-1);

%% Old code
%     N = numel(v);
% 
%     v = sort(v(:));
%     
%     x = zeros(N,1);
%     y = zeros(N,1);
%     x(1) = v(1);
%     y(1) = 1;
%     ind = 1;
%     for ii = 1:(N-1)
%         if v(ii+1)>v(ii)+ eps(10);
%             x(ind+1) = v(ii+1);
%             y(ind+1) = y(ind)+1;                
%             ind = ind + 1;
%         else
%             y(ind) = y(ind) + 1;
%             if ii == N-1
%                 x(ind) = v(ii);
%             end
%         end
%         
%     end
%     
%     [~,~,x] = find(x);
%     [~,~,y] = find(y);
%     
%     y = y/N;
%     



end