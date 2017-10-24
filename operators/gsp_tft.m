function Xdot = gsp_tft(G,X)
% GSP_TFT Compute the time Fourier Transform of a time-vertex signal
%   Usage:  X = gsp_tft(G,Xhat)
%
%   Input parameters:
%         G      : Time-Vertex graph structure
%         X      : Time-Vertex signal
%   Output parameters:
%         Xdot   : Time-Vertex Time Fourier Coefficients
%
%   Compute the time Fourier Transform of a time-vertex signal

% Author : Francesco Grassi
% Date   : September 2016


NFFT = G.jtv.NFFT;

if isempty(NFFT)
    NFFT = size(X,2);
end;
    

normalize = 1/sqrt(NFFT);

switch G.jtv.transform
    case 'dft'
        Xdot = fft(X,NFFT,2)*normalize;
    case 'dct'
        if size(X,3)>1
            Xdot = ipermute(dct(permute(X,[2 1 3]),NFFT),[2 1 3]);
        end
        Xdot = dct(X.',NFFT).';
    
    otherwise
        error('Unknown transform');
end
