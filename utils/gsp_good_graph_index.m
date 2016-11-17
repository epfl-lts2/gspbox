function [ index ] = gsp_good_graph_index( G, X, param )
%GSP_GOOD_GRAPH_INDEX Index testing how well a given graph G, matches some data X
%   Usage: gsp_good_graph(G, X);
% 
%   Input parameters:
%       G          : the graph
%       X          : a data matrix
%       param      : structure of optional parameters
%   Output parameters: 
%       index      : the computed index 
% 
%   A wrapper function with which one may test how well a given graph G,
%   matches some data X.
% 
%   Example:
%       G = gsp_2dgrid(16);
%       X = pinv(full(G.L))  randn(G.N, G.N);
%       param.verbose = 1;
%       param.index = 'tcer';
%       index = gsp_good_graph_index(G, X, param)
%       param.index = 'stationarity';
%       index = gsp_good_graph_index(G, X, param)
%
%   Optional paramaters
%   -------------------
%
%    param.index*: 'tcer' or 'stationarity' (default 'tcer'). 
% 
%   See also: gsp_learn_tcer, gsp_stationarity_ratio
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_good_graph_index.php

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


% Author  : Andreas Loukas
% Date    : 15 Nov 2016

% Handle input
if nargin < 3, param = struct(); end
if not(isfield(param, 'index'));   param.index = 'tcer'; end;

switch param.index,
  
    case 'tcer',    
        index = gsp_learn_tcer(G, X, param);
        
    case 'stationarity'        
        index = gsp_stationarity_ratio(G, X*X', param);
        
    otherwise, 
        error('uknown index.');

end
end


