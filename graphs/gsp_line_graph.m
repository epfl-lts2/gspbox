function LG = gsp_line_graph(G,param)
%GSP_LINE_GRAPH Create Line Graph (or edge-to-vertex dual graph) of graph G
%   Usage: LG = gsp_line_graph(G);
%          LG = gsp_line_graph(G,param);
%
%   Input parameters:
%       G        : Graph structure
%       param    : Structure of optional parameter
%   Output parameters:
%       LG       : Graph structure
%
%  
%   Use param.coords  to compute coordinates using halfway point of
%   coordinates of graph G (if available)
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/graphs/gsp_line_graph.html

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
 
% Author: Francesco Grassi
% Date  : July 2016
 
 
if nargin<2
    param=struct;
end
 
if ~isfield(G,'B'); G = gsp_incidence(G); end;
if ~isfield(param,'coords'); param.coords=isfield(G,'coords'); end;
 
 
%Line Graph adjacency matrix
W_LG = G.B*G.B';
W_LG = W_LG - diag(diag(W_LG));
 
LG = gsp_graph(W_LG);
 
 
if param.coords
    % coordinates
    halfway = @(x) squareform(bsxfun(@plus,x,x')-diag(2*x))/2;
    ndim = size(G.coords,2);
    coords = zeros(G.Ne,ndim);
    
    %small loop over the graph dimension (2 or 3 iterations)
    for n=1:ndim
        diff_coord=halfway(G.coords(:,n));
        diff_coord = diff_coord(squareform(G.A)>0);
        coords(:,n) = diff_coord;
    end
    
    LG.coords = coords;
    
end
 
 
 
 
 



