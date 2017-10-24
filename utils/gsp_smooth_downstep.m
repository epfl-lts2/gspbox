function y = gsp_smooth_downstep(x, a, o)
%GSP_SMOOTH_DOWNSTEP Smooth downstep function from 1 to 0
%   Usage: y = gsp_smooth_downstep(x);
%          y = gsp_smooth_downstep(x, a);
%          y = gsp_smooth_downstep(x, a, o);
%   
%   Input parameters:
%          x        : input value
%          a        : smoothing parameter (default 1)
%          o        : offset (default 1)
%   Ouput parameters:
%          y        : output value of the function
%
%   This function is a smooth downstep from 1 to 0 arround *o* using the
%   following function:
%   
%   .. f(x) = s(o-x)
%
%   .. math:: f(x) = s(o-x)
%   
%   where 
%
%   ..          /   0                                      if x < -1
%   ..  s(x) = | exp(-a/x) / ( exp(-a/x) + exp(-a/(1-x)) ) if x in [-1, 1]
%   ..          \   1                                      if x > 1
%
%   .. math:: s(x)=\begin{cases} 0 & \mbox{if }x<-1 \\ \frac{e^{-\frac{a}{x}}}{e^{-\frac{a}{x}}+e^{-\frac{a}{1-x}}} & \mbox{if }x\in[-1,1]\\ 1 & \mbox{if }x>1 \end{cases}
%

% Author: Nathanael Perraudin
% Date  : 5 mai 2016


if nargin<3
    o = 1;
end
    
if nargin<3
    a = 1;
end

y = gsp_smooth_step(o-x, a);

end


