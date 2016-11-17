function L = gsp_laplacian(W, lap_type)
%GSP_LAPLACIAN compute the graph Laplacian from undirected adjacency matrix
%   Usage: L = gsp_laplacian(W, laplacian_type);
%          G = gsp_laplacian(W);
%
%   Input parameters:
%       W   : Weighted adjacency matrix
%       laplacian_type: Type of laplacian: 'combinatorial' or 'normalized'
%   Output parameters:
%       L   : Graph Laplacian
%
%   This function creates the graph laplacian given a weighted adjacency
%   matrix.
%
%   The variable laplacian_type contains the different laplacian types.
%   Available laplacian types:
%
%    combinatorial*: Non normalized laplacian (default).
%
%          L =  D  - W 
%
%   For directed graphs, see gsp_create_laplacian.
%
%   References:
%     F. Chung. Laplacians and the cheeger inequality for directed graphs.
%     Annals of Combinatorics, 9(1):1--19, 2005.
%     
%     
%
% see also: gsp_create_laplacian
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_laplacian.php

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

% Author: Vassilis Kalofolias
% Date  : June 2016


if nargin < 2
    lap_type = 'combinatorial';
end

n = size(W,1);

d = sum(W,2);
switch lap_type
    case 'combinatorial'
        L = diag(d)-W;
    case 'normalized'
        if issparse(W)
            Dn = diag(d.^(-0.5));
            % CAREFUL: a zero row/column in W is followed by an inf in d,
            % but their multiplication gives zero again (because of using
            % sparse).
            L = speye(n) - Dn * W * Dn;
        else
            ind = d>0;
            L = - W;
            % Ln = D^(-1/2) L D^(-1/2)
            L(ind, :) = bsxfun(@times, L(ind, :), 1./sqrt(d(ind)));
            L(:, ind) = bsxfun(@times, L(:, ind), 1./sqrt(d(ind))');
            % put back diagonal to identity
            % Note: for disconnected nodes we should still have 1 on diagonal
            % (limit of L for W -> 0)
            L(1:n+1:end) = 1;
    %         L(ind, ind) = speye(nnz(ind))-Dn*W(ind, ind)*Dn;
        end
    case 'randomwalk'
        L = speye(n) - diag(d.^(-1))*W;

    case 'none'
        L = sparse(0);
    otherwise
        error(' Unknown laplacian type')
end

