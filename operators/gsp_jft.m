function Xhat = gsp_jft(G,X)
%GSP_JFT Compute the Joint time-vertex Fourier Transform
%   Usage:  Xhat = gsp_jft(G,X)
%
%   Input parameters:
%         G          : Time-Vertex graph structure
%         X          : Time-Vertex signal
%   Output parameters:
%         Xhat       : Joint Time-Vertex Fourier Transform of X
%
%
%   'gsp_jft(G,X)' computes the joint time-vertex Fourier transform of the time-vertex
%    signal X with respect to the Fourier basis of the graph G: U_G and the DFT basis U_T.
%
%      X_hat = U_G' * X * conj(U_T)
%
%
%   To compute the Fourier basis of a graph G, you can use the function:
%
%           G = gsp_compute_fourier_basis(G);
%
%   Example:
%
%           N = 50; T=200;
%           G = gsp_sensor(N);
%           G = gsp_jtv_graph(G,T);
%           G = gsp_compute_fourier_basis(G);
%           X = sin((1:N)'*(1:T)*pi/(4*N));
%           Xhat = gsp_jft(G,X);
%           imagesc(fftshift(abs(Xhat),2));
%
%   See also: gsp_ijft, gsp_compute_fourier_basis
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/operators/gsp_jft.html

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

% Author : Francesco Grassi
% Date   : September 2016


if isempty(G.jtv.NFFT)
    NFFT = size(X,2);
else
    NFFT = G.jtv.NFFT;
end

normalize = 1/sqrt(NFFT);

switch G.jtv.transform
    case 'dft'
        Xhat = fft(gsp_gft(G,X),NFFT,2)*normalize;
    case 'dct'
        if size(X,3)>1
            Xhat = ipermute(dct(permute(gsp_gft(G,X),[2 1 3]),NFFT),[2 1 3]);
        else
            Xhat = dct(gsp_gft(G,X).',NFFT).';
        end
        
    otherwise
        error('Unknown transform');
end

end
