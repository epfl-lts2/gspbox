function [y,x] = gsp_ddf2dcdf(v,x)
%GSP_DDF2DCDF Discrete density function to dicrete cumulative density function
%   Usage : [y] = gsp_ddf2dcdf(v);
%           [y,x] = gsp_ddf2dcdf(v,x);
%
%   Input parameters:
%       v   : Discrete density function
%       x   : Axis value (Optional)
%   Output parameters:
%       y   : Discrete cumulative ensity function
%       x   : Axis value (Same as input)
%
%   This function goes from discrete density function to discrete
%   cumulative density function.
%   


N = length(v);
y = zeros(N,1);
y(1) = v(1);
for ii = 1:(N-1) 
    y(ii+1) = y(ii) +v(ii+1);
end

y = y/y(end);

end