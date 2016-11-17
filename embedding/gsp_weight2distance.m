function D = gsp_weight2distance(G, method, sigma)
%GSP_WEIGHT2DISTANCE Distance matrix from weight matrix
%   Usage: D = gsp_weight2distance(G);
%          D = gsp_weight2distance(G, method);
%          D = gsp_weight2distance(G, method, nnn);
%
%   Input parameters
%         G         : Graph or weight matrix
%         method    : Type of kernel used to compute the weight matrix
%         sigma     : kernel paramter
%
%
%   Output parameters
%         D         : Distance matrix (sparse)
%
%   This function uses the weight matrix of a graph in order to compute the
%   distance matrix D. If G is a structure it is considered as a graph
%   otherwise it is considered as a weight matrix. The resulting matrix is
%   in sparse form. The indices idx are of the nearest neighbors or (if k
%   is not used) the indices of all non zero elements of the weight matrix.
%   In both cases indices are saved in a cell where idx{ii} contains a
%   vector jj of indices where D(ii, jj) are the Distances of the nearest
%   neighbours of element ii.
%
%   The type of method one can use are the following:
%    'exp'        : exponential kernel (e^(frac{-d_{ij}}{sigma^2}))
%    '1/x'        : inverse of x kernel (frac{1}{sigma+d_{ij}})
%    '1/x^2'      : inverse of x^2 kernel (frac{1}{(sigma+d_{ij})^2})
%    'resistance' : Resistance distance.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/embedding/gsp_weight2distance.php

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

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015



if isstruct(G)
    W = G.W;
else
    W = G;
end

if nargin<2
    method = 'exp';
end
if nargin<3
    sigma = 1;
end

switch lower (method)
    case 'exp'
        if sum(W(:)>1)
            error('W cannot have entries bigger than 1')
        end
        %             d =  -  log(w) * sigma;
        %             D = sparse(rows,cols,d);
        D = spfun(@(w) sqrt(abs(- log(w))) * sigma, W);
        % w = 1/d  =>
    case '1/x'
        if any(W(:) >= sigma)
            error('W cannot have entries bigger than 1/sigma')
        end
        %             d = 1./w - sigma;
        %             D = sparse(rows,cols,d);
        D = spfun(@(w) 1./w - sigma, W);
    case '1/x^2'
        if any(W(:) >= sigma^2)
            error('W cannot have entries bigger than 1/sigma^2')
        end
        %             d = sqrt(1./w) - sigma;
        %             D = sparse(rows,cols,d);
        D = spfun(@(w) sqrt(1./w) - sigma, W);
        
    case 'resistance'
        D = gsp_resistance_distances(G);
        %     d = d(rows,cols);
        D = sparse(D);
    otherwise
end


end
