function jtwgft = gsp_jtwgft(G, w, s, param)
%GSP_JTWGFT Compute the joint vertex-time windowed graph fourier transform
%   Usage:  jtwgft = gsp_jtwgft(G, w, s)
%           jtwgft = gsp_jtwgft(G, w, s, param)
%
%   Input parameters:
%         G          : Graph structure
%         w          : windows
%         s          : vertex-time signal
%         param      : structure of optional parameters
%   Output parameters:
%         jtwgft     : Joint vertex-time windowed graph Fourier transform
%
%   Additional parameters
%   ---------------------
%    param.boundary    : 'periodic' or 'reflecting' (default 'periodic')
%    param.a           : hop size in time (default 1)
%    param.M           : number of frequencies (default size(s,2))
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_jtwgft.php

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

% Author : Nathanael Perraudin
% Date   : 30 April 2016

if nargin<4
    param=struct;
end
if ~isfield(param,'boundary'), param.boundary = 'periodic'; end
if ~isfield(param,'a'), param.a = 1; end
if ~isfield(param,'M'), param.M = size(s,2); end


switch param.boundary
    case 'periodic'
        % 1) perform DGT
        S = dgt(transpose(s),w,param.a,param.M);
        % 2) Reshape
        S = reshape(S,[],G.N);
        % 3) GFT
        Shat = transpose(gsp_gft(G,transpose(S)));
        % 4) Reshape
        jtwgft = reshape(Shat,param.M,[],G.N);
    case 'reflecting'
        % 1) perform DGT
        S = wmdct(transpose(s),w,param.M,size(s,2));
        % 2) Reshape
        S = reshape(S,[],G.N);
        % 3) GFT
        Shat = transpose(gsp_gft(G,transpose(S)));
        % 4) Reshape
        jtwgft = reshape(Shat,param.M,[],G.N);
    otherwise
        error('Unknown boundary condition');
end

end
