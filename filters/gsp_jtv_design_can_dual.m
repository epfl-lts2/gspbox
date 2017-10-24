function [ gd,filtertype ] = gsp_jtv_design_can_dual(g,filtertype)
%GSP_JTV_DESIGN_CAN_DUAL This function returns the canonical dual of the time-vertex filterbank g
%   Usage:  [gd,filtertype] = gsp_jtv_design_can_dual( g,filtertype );
%
%   Inputs parameters:
%       g          : cell array of time-vertex filters
%       filtertype : Filter domain (ts,js,ts-array,js-array)
%
%   Ouputs parameters:
%       g          : cell array of canonical dual time-vertex filters
%       filtertype : Filter domain (ts,js)
%
%   This function returns the canonical dual of the time-vertex filterbank g
%

% Author: Francesco Grassi
% Date:   September 2016

Nf = size(g,1);
gd = cell(Nf,1);

for n = 1:Nf
    gd{n} = @(x,t) can_dual(g,filtertype,n,x,t);
end

switch filtertype
    case {'ts','ts-array'}
        filtertype = 'ts';
    case {'js','js-array'}
        filtertype = 'js';
end

end


function sol = can_dual(g,ft,n,x,t)


if ~isvector(x);x=x(:,1);end
if ~isvector(t);t=t(1,:);end

sol = gsp_jtv_evaluate_can_dual( g,ft,x,t,n );


end

