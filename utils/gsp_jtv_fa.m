function [f, label] = gsp_jtv_fa(G,shift)
%GSP_JTV_FA Frequency axis for the joint time vertex framework
%   Usage:  f = gsp_jtv_fa(G);
%
%   Input parameters:
%       G       : Time-vertex graph structure
%       shift   : Boolean value: 1 to apply fftshift on the frequency vector, 0 otherwise
%   Ouput parameters:
%       f       : Frequency axis (row vector)
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_jtv_fa.php

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

% Author: Nathanael Perraudin, Francesco Grassi
% Date  : September 2016

if nargin<2
    shift = 0;
end

if ~gsp_check_jtv(G)
    error('GSP_JTV_FA needs the time dimension. Use GSP_JTV_GRAPH');
end

if isempty(G.jtv.NFFT)
    NFFT = G.jtv.T;
else
    NFFT = G.jtv.NFFT;
end

if shift
    f = fftshift(gsp_cfa( NFFT,G.jtv.fs )');
else
    f = gsp_cfa( NFFT,G.jtv.fs )';
end




label = '\omega';

end
