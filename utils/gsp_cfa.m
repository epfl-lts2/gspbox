function [ fa ] = gsp_cfa( T,fs )
%GSP_CFA Create frequency axis
%   Usage: fa = gsp_cfa(N);
%          fa = gsp_cfa(N,fs);
%
%   Input parameters:
%       N   : Number of samples
%       fs  : Sampling frequency (default 1)
%   Ouput parameters:
%       fa  : Frequency axis
%   
%   This function create a normalized frequency axis corresponding to the
%   frequency of the function fft.
%

% Author: Nathanael Perraudin

if nargin<2
    fs = 1;
end

if mod(T,2)
    fa = linspace(-0.5+1/(2*T),0.5-1/(2*T),T);
else
    fa = linspace(-0.5,0.5-1/T,T);
end

fa = fa*fs;

fa = ifftshift(fa);
fa = fa(:);

end
