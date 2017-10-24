function [sol] = gsp_wiener_inpainting_exact(G,y, M, psd, psd_noise)
%GSP_WIENER_INPAINTING_MATRIX Exact solution for wiener in-painting
%   Usage:  sol = gsp_wiener_inpainting_exact(G,y, M, psd, psd_noise)
%           sol = gsp_wiener_inpainting_exact(G,y, M, psd, psd_noise, param)
%
%   Input parameters:
%         G          : Graph (GSP structure)
%         y          : Measurements (column vector)
%         M          : Mask (vector)
%         psd        : PSD filter (anonymous function)
%         psd_noise  : PSD filter of the noise or single number
%   Output parameters:
%         sol        : Solution
%         infos      : Convergence informations
%
%   This function returns the exact solution to:
%
%        argmin_x || M x - y ||_2^2 + || w(L) x ||_2^2 
%
%   Please refer to the reference for more information about this problem.
%
%   References:
%     N. Perraudin and P. Vandergheynst. Stationary signal processing on
%     graphs. arXiv preprint arXiv:1601.02522, 2016.
%     
%     
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/stationarity/gsp_wiener_inpainting_exact.html

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

% Author : Nathanael Perraudin
% Date: 6 January 2016


if ~isnumeric(psd_noise) || numel(psd_noise)>1
    error('The psd of the noise has to be a scalar');
end

if ~gsp_check_fourier(G)
    error('You need the Fourier basis to use this function.');
end

psdL = gsp_filter_analysis(G,psd,eye(G.N));

if (numel(M) == size(M,1)) || (numel(M) == size(M,2))
    indl = find(M);
    indu = find(1-M);        
else   
    error('I cannot handle this case yet');
end

if psd_noise==0
    sol = zeros(G.N,size(y,2));
    sol(indl,:) = y(indl,:); 
    sol(indu,:) = psdL(indu,indl) * ...
        (( psdL(indl,indl) + eye(numel(indl)) * eps(10*max(psd(G.e)) )) \ y(indl,:));
else
%     Mop =@(x) bsxfun(@times,M,x);
%     sol = ( psdL * diag(M) + eye(G.N) * psd_noise ) \ (psdL * Mop(y) );
    sol = psdL(:,indl) * ...
        ( (psdL(indl,indl)+ eye(numel(indl))* psd_noise ) \  y(indl,:));
end

    

end
