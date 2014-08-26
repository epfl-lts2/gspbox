function rd = gsp_resistance_distance(G,param)
%GSP_RESISTANCE_DISTANCE : Compute the resitance distances of a graph
%   Usage: rd = gsp_resistance_distance(G);
%          rd = gsp_resistance_distance(L);
%
%   Input parameters:
%       G    : Graph structure or Laplacian matrix (L)
%       param: optional parameters
%   Output parameters:
%       rd   : distance matrix
%
%   This function compute the resistance distance of a graph. The distance
%   between two nodes is defined as the inverse of the weight matrix. For
%   example the distance matrix:
%
%           dist = [0, 3, 1;...
%                   3, 0, 2;...
%                   1, 2, 0];
%
%   Correspond to the weight matrix:
%
%           W = [0, 1/3, 1/1;...
%                1/3, 0, 1/2;...
%                1/1, 1/2, 0];
%
%   The function will compute the resistance distance following the
%   Kirshoff's law. In the our example it is:
%
%           rd2 = [0, 3/2, 5/6;...
%                  3/2, 0, 4/3;...
%                  5/6, 4/3, 0]
%
%   In matlab, you can reprocude this example using:
%
%           % The weigh 
%           dist = [0, 3, 1;...
%                   3, 0, 2;...
%                   1, 2, 0];
%           % The weight is the inverse of the distance...
%           W = dist.^(-1);
%           % Fix the diagonal
%           W([1,5,9])=0;    
%           G = gsp_graph(W);
%           rd = gsp_resistance_distance(G)
%           % Resitance computed by hand
%           rd2 = [0, 3/2, 5/6;...
%                  3/2, 0, 4/3;...
%                  5/6, 4/3, 0]
%
%   param is an optional structure that contains the following field
%
%    param.verbose*: display parameter - 0 no log - 1 display warnings
%     (default 1)
%   
%   References:
%     D. J. Klein and M. Randić. Resistance distance. Journal of Mathematical
%     Chemistry, 12(1):81-95, 1993.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_resistance_distance.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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

% Author: Nathanael Perraudin

if nargin < 2
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1;

if isstruct(G)
    % Use the non normalized laplacian
    if ~strcmp(G.lap_type, 'combinatorial')
        G = gsp_create_laplacian(G,'combinatorial');
        if param.verbose
            fprintf(['Compute the combinatorial laplacian ',...
            'for the resitance distance\n']);
        end
    end
    L = G.L;
else
    L = G;
end



pseudo=pinv(full(L));
N = size(L,1);

d = diag(pseudo);

rd = repmat(d,1,N) + repmat(d',N,1) - pseudo - pseudo';

end



