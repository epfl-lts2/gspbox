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
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/operators/gsp_tft.html

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

