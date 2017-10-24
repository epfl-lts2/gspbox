function [ h ] = gsp_jtv_evaluate_can_dual( g,filtertype,x,t,n )
%GSP_JTV_EVALUATE_CAN_DUAL Evaluates the canonical dual of time-vertex filterbank
%   Usage: h = gsp_jtv_evaluate_can_dual( g,x,t )
%          h = gsp_jtv_evaluate_can_dual( g,x,t,n )
%
%   Inputs parameters:
%       g       : Cell array of time-vertex filters
%       filtertype : Filter domain (ts,js,ts-array,js-array)
%       x       : Data
%       t       : Time
%       n       : Index of filter in the filterbank to evaluate (default 0 = all filters)
%
%   Ouputs parameters:
%       h       : Cell array of dual time-vertex filters
%
%   Evaluates the canonical dual of time-vertex filterbank
%

% Author: Francesco Grassi
% Date: July 2016

if ~gsp_check_filtertype(filtertype)
    error('Invalid filtertype');
end

if nargin<5
    n = 0;
end

N = length(x);
T = length(t);
Nf = size(g,1);

% Compute coefficient of g
gcoeff = gsp_jtv_filter_evaluate(g,filtertype,x,t);
gcoeff = reshape(gcoeff,N*T,Nf);

if n
    h = sum(conj(gcoeff).*gcoeff,2).^-1.*gcoeff(:,n);
else
    h = zeros(N*T,Nf);
    for ii=1:Nf
        h(:,ii)= sum(gcoeff.^2,2).^-1.*gcoeff(:,ii);
    end 
end

if any(strcmpi(filtertype,{'ts','ts-array'}))
    h = ifft(reshape(h,N,T,[]),[],2)*sqrt(T);
else
    h = reshape(h,N,T,[]);
end

