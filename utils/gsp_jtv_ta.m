function [t,label] = gsp_jtv_ta(G,lag)
%GSP_JTV_TA Time axis for the joint time vertex framework
%   Usage:  t = gsp_jtv_ta(G);
%
%   Input parameters:
%       G       : Time-vertex graph structure
%       lag     : Compute lag axis
%   Ouput parameters:
%       t       : Time or Lag axis (row vector)
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_jtv_ta.php

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

% Author: Francesco Grassi
% Date  : September 2016

if nargin<2
    lag = 0;
end

if isfield(G.jtv,'T')
    T = G.jtv.T;
    fs = G.jtv.fs;
    if lag
        t = -(T-1)/fs:1/fs:(T-1)/fs;
        label = 'Lag';
    else
        t = 0:1/fs:(T-1)/fs;
        label = 'Time';
    end
else
    error('GSP_JTV_TA needs the time dimension. Use GSP_JTV_GRAPH');
end

end
