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
%     .. argmin_x || M x - y ||_2^2 + || w(L) x ||_2^2 
%
%     .. math:: arg\min_x \| M x - y \|_2^2 + \| w(L) x \|_2^2 
%
%   Please refer to the reference for more information about this problem.
%
%   References: perraudin2016stationary

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