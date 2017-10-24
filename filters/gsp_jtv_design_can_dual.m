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
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_design_can_dual.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

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


