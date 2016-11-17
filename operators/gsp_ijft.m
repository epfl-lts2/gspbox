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
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_ijft.php

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

