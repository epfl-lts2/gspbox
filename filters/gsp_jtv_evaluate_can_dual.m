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
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_evaluate_can_dual.html

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


