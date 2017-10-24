function X = gsp_ijft(G,Xhat)
% GSP_IJFT Compute the inverse Joint time-vertex Fourier Transform
%   Usage:  X = gsp_ijft(G,Xhat)
%
%   Input parameters:
%         G          : Time-Vertex graph structure
%         Xhat       : Time-Vertex Fourier coefficients
%   Output parameters:
%         X          : Time-Vertex signal
%
%   Compute the inverse Joint time-vertex Fourier Transform

% Author : Francesco Grassi
% Date   : September 2016


if isempty(G.jtv.NFFT)
    NFFT = size(Xhat,2);
else
    NFFT = G.jtv.NFFT;
end

normalize = sqrt(NFFT);

switch G.jtv.transform
    case 'dft'
        X = ifft(gsp_igft(G,Xhat),NFFT,2)*normalize;
    case 'dct'
        if size(Xhat,3)>1
            error('Not implemented')
        end
        X = idct(gsp_igft(G,Xhat).',NFFT).';
    otherwise
        error('Unknown transform');
end


end
