function varargout = histnorm(varargin)
% HISTNORM Histogram normalized
%   [...] = HISTNORM(...) works like HIST, but the frequency is normalized
%   so that area sum is 1.
%
%   Bonus usage!
%   [...] = HISTNORM(..., 'plot') plots and returns the output arguments. 
%   Be sure 'plot' is the last argument.
%
%   Example:
%       data = randn (10000, 1);
%       [xo,no] = histnorm(data, 101, 'plot');
%       hold on
%       plot (no, normpdf(no), 'r');
%       hold off
%
%   See also: HIST.
%
%   Written by Arturo Serrano and downloaded from MATLAB Central File
%   Exchange
%
%   Copyright 2009 DWTFYW.
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/utils/histnorm.php

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

% check whether plot is done
doPlot = 0;
if ischar (varargin{end}) && strcmpi ('plot', varargin{end})
    doPlot = 1;
    varargin = varargin(1:end-1);
elseif nargout == 0
    doPlot = 1;    
end

% normalize so the "area" is 1
[xo,no] = hist (varargin{:});
binwidths = diff ([min(varargin{1}) no(1:end-1)+diff(no)/2 max(varargin{1})]);
xonorm = xo/sum (xo .* binwidths);
varargout = {xonorm, no};
varargout = varargout(1:nargout);
figure
hold on
% do plot
if doPlot
    cax = axescheck(varargin{:});

%     % bored way: bar plot
%     if isempty (cax)
%         bar (no, xonorm, 'hist');
%     else
%         bar (cax, no, xonorm, 'hist');
%     end
    
    % funny way: modify vertices of bar plot
    hist (varargin{:});
    if isempty (cax)
        cax = gca;
    end
    ch = findobj (get (cax, 'children'), 'type', 'patch'); ch = ch(1);
    vertices = get (ch, 'vertices');
    for idx = 1:numel (xonorm)
        vertices((idx-1)*5+[3 4],2) = xonorm(idx);     % hope it works :)
    end
    set (ch, 'vertices', vertices);
end

