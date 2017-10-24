function y = gsp_smooth_step(x, a)
%GSP_SMOOTH_STEP Smooth step function from 0 to 1 arround 0.5
%   Usage: y = gsp_smooth_step(x);
%          y = gsp_smooth_step(x, a)
%   
%   Input parameters:
%          x        : input value
%          a        : smoothing parameter (default 1)
%   Ouput parameters:
%          y        : output value of the function
%
%   This function is a smooth step from 0 to 1 arround 0 using the
%   following function:
%
%   ..          /   0                                      if x < 0
%   ..  s(x) = | exp(-a/x) / ( exp(-a/x) + exp(-a/(1-x)) ) if x in [0, 1]
%   ..          \   1                                      if x > 1
%
%   .. math:: s(x)=\begin{cases} 0 & \mbox{if }x<0 \\ \frac{e^{-\frac{a}{x}}}{e^{-\frac{a}{x}}+e^{-\frac{a}{1-x}}} & \mbox{if }x\in[0,1]\\ 1 & \mbox{if }x>1 \end{cases}
%


% Author: Nathanael Perraudin
% Date  : 5 mai 2016


if nargin<2
    a = 1;
end

y=fx(x,a);
y = y./(y+fx(1-x,a));

end

function y=fx(x,a)
    y = exp(-a./x);
    y(x<0)=0;
end
