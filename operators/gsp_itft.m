function X = gsp_itft(G,Xdot)
% GSP_ITFT Compute the inverse Time Fourier Transform of a time-vertex signal
%   Usage:  X = gsp_itft(G,Xhat)
%
%   Input parameters:
%         G      : Time-Vertex graph structure
%         Xdot   : Time-Vertex Time Fourier Coefficients
%   Output parameters:
%         X      : Time-Vertex signal
%   Compute the inverse Time Fourier Transform of a time-vertex signal

% Author : Francesco Grassi
% Date   : September 2016


NFFT = G.jtv.NFFT;

if isempty(NFFT)
    NFFT = size(Xdot,2);
end;
    

normalize = sqrt(NFFT);

switch G.jtv.transform
    case 'dft'
        X = ifft(Xdot,NFFT,2)*normalize;
        if sum(abs(vec(imag(X)))) < 1e-13 *norm(X(:))
            X = real(X);
        end
    case 'dct'
        X = idct(Xdot.',NFFT).';
    otherwise
        error('Unknown transform');
end


end
