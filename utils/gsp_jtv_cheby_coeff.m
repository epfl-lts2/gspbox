function c = gsp_jtv_cheby_coeff(G, filter, filtertype, m, N,param)
%GSP_JTV_CHEBY_COEFF : Compute Chebyshev coefficients for a time-vertex filterbank
%   Usage: c = gsp_jtv_cheby_coeff(G, filter);
%          c = gsp_jtv_cheby_coeff(G, filter, m);
%          c = gsp_jtv_cheby_coeff(G, filter, m, N);
%
%   Input parameters:
%       G       : Time-Graph structure
%       filter  : Cell array of time-vertex filters
%       filtertype : Filter domain (ts,js,ts-array,js-array)
%       m       : maximum order Chebyshev coefficient to compute (default 30)
%       N       : grid order used to compute quadrature (default is m+1)
%       param   : structure of optional parameter
%   Output parameters
%       c       : matrix of Chebyshev coefficients
%
%   Additional parameters
%   --------------------- 
%    param.verbose : Verbosity level (1 display the warning - 0 no log)
%     (default 1).
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_jtv_cheby_coeff.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

%
% Author: Francesco Grassi
% Date   : July 2016

if nargin < 6
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;

if nargin < 4
    m = 30;
end

if nargin < 5
    N = m+1;
end

if isstruct(G)
    if or(~isfield(G.jtv,'T'),~isfield(G.jtv,'fs'));
        error('GSP_JTV_CHEBY_COEFF need time dimension. Use GSP_JTV_GRAPH.')
    end
    
else
    error('Invalid graph structure.')
end


T = G.jtv.T;

if G.jtv.extension
    tau = 2*G.jtv.T-1;
else
    tau = T;
end


if iscell(filter)
    Nf = length(filter);
    c = zeros(m+1,tau,Nf);
    for ii = 1: Nf
        c(:,:,ii) = gsp_jtv_cheby_coeff(G, filter{ii},filtertype, m, N,param);
    end
    return;
end

if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);
    if param.verbose
        warning(['GSP_JTV_CHEBY_COEFF: The variable lmax is not ',...
            'available. The function will compute it for you. ',...
            'However, if you apply many time this function, you ',...
            'should precompute it using the function: ',...
            'gsp_estimate_lmax']);
    end
end


arange = [0, G.lmax];



a1=(arange(2)-arange(1))/2;
a2=(arange(2)+arange(1))/2;



t = gsp_jtv_ta(G);


c = zeros(N,T);

param.domain = 'joint-spectral';
for ii=1:m+1
    c(ii,:) = sum( gsp_jtv_filter_evaluate(filter,filtertype,a1* cos( (pi*((1:N)-0.5))/N) + a2,t,param ) .* ...
        repmat(cos( pi*(ii-1)*((1:N)-0.5)/N).',1,T)) *2/N;
end

if G.jtv.extension
    c = [c zeros(N,T-1)];
end

end








