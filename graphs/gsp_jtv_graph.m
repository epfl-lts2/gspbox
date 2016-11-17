function G = gsp_jtv_graph(G,T,fs,param)
% GSP_JTV_GRAPH Add time information to the graph structure
%   Usage:  G = gsp_jtv_graph(G);
%           G = gsp_jtv_graph(G,T);
%           G = gsp_jtv_graph(G,T,fs);
%           G = gsp_jtv_graph(G,T,fs,param);
%
%   Input parameters:
%         G          : Graph structure
%         T          : Time length
%         fs         : Sampling frequency (default 1)
%         param      : Structure of optional parameters
%
%   Output parameters:
%         G          : Time-Vertex Graph structure
%
%   This function adds the time domain information to the structure G.
%   The fields stored inside G.jtv are:
%   - T         : Time length
%   - fs        : Sampling frequency
%   - NFFT      : Number of frequency point (default [])
%   - Transform : Time Fourier basis: 'dft' or 'dct'
%   - Extension : Signal will be zero padded in time to take into account negative lag
%   - Lag       : Lag axis
%   - omega     : Frequency axis
%   - DiffT     : Time Gradient
%   - LT        : Time Laplacian
%
%   Additional parameters
%   ---------------------
%
%    param.transform  : Fourier basis: 'dft' or 'dct'. (default 'dft')
%    param.approx     : Stencil approx for the gradient: 'forward' or 'backward'. (default 'forward')
%    param.extension  : Signal will be zero padded in time when needed. (default '0')
%    param.NFFT       : Number of frequency point. (default [])
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_jtv_graph.php

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

% Author :  Francesco Grassi
% Date   : September 2016

if nargin<4
    param = struct;
end

if ~isstruct(G)
    error('G is not a valid graph');
end

if nargin<2 || ~isnumeric(T) || T<1
    error('Time length T must be a numeric value greater than 0');
end

if nargin<3 || isempty(fs)
    fs=1;
end

if ~isnumeric(fs) || fs<0
    error('Sampling frequency fs must be a numeric value greater than 0');
end

if ~isfield(param,'NFFT'),      param.NFFT = []; end
if ~isfield(param,'transform'), param.transform = 'dft'; end
if ~isfield(param,'approx'),    param.approx = 'forward'; end
if ~isfield(param,'extension'), param.extension = 0; end


%% JTV PARAMETERS
G.jtv = struct;
G.jtv.T = T;
G.jtv.fs = fs;
G.jtv.NFFT = param.NFFT;
G.jtv.transform = param.transform;
G.jtv.extension = param.extension;
G.jtv.omega = gsp_jtv_fa(G);

if G.jtv.extension
    G.jtv.lag = 2*G.jtv.T-1;
else
    G.jtv.lag = T;
end

%% DIFFERENTIAL OPERATORS
z = ones(T,1);

G.jtv.LT = spdiags([-z 2*z -z], [-1 0 1], T, T);

G.jtv.DiffT = spdiags([-z z], [0 1], T, T);

switch param.transform
    case 'dft'
        G.jtv.DiffT(T,1) = 1;
    case 'dct'
        G.jtv.DiffT(T,T-1) = 1;
    otherwise
        error('Unknown transform');
end
 
if strcmpi(param.approx,'backward')
    G.jtv.DiffT = -G.jtv.DiffT.';
end

% toeplitz([-1 zeros(T-2,1) 1],[-1 1 zeros(T-2,1)]); %periodic forward
% toeplitz([1 -1 zeros(T-2,1)],[1 zeros(T-2,1) -1]); %periodic backward (= -transpose(periodic forward))


end

