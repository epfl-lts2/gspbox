function coords = gsp_laplacian_eigenmaps(G, dim, param)
%GSP_LAPLACIAN_EIGENMAPS Laplacian eigenmaps embedding 
%   Usage: coords = gsp_laplacian_eigenmaps(G, dim);
%          coords = gsp_laplacian_eigenmaps(G, dim, param);
%
%   Input parameters
%         G      : Graph
%         dim    : Dimension of the embedding
%         param  : Structure of optional paramters
%   Output parameters
%         coords : Coordinates of the embedding
% 
%   This function uses the weight matrix of a graph G, in order to compute
%   a dim -dimensional embedding (output coordinates). The algorithm used
%   is Laplacian eigenmaps. Warning, this function might not work if the
%   graph is not connected.
%
%   param is a structure optional parameters:
%     
%    param.tol : Tolerance for the spectral gap (default 1e-6).   
%
%
%   See also: gsp_update_coordinates gsp_isomap gsp_lle
%   
%   Demo: gsp_demo_graph_embedding
%
%   References:
%     M. Belkin and P. Niyogi. Laplacian eigenmaps and spectral techniques
%     for embedding and clustering. In NIPS, volume 14, pages 585--591, 2001.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/embedding/gsp_laplacian_eigenmaps.php

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


if nargin<3
    param = struct;
end

if ~isfield(param,'tol'), param.tol = 1e-6; end


% Laplacian Eigenmaps
options.disp = 0;
options.isreal = 1;
options.issym = 1;
% only need bottom (no_dims + 1) eigenvectors
if strcmp(G.type,'normalized') % Bug G.type
    [coords, e] = eigs(G.L, dim + 1, param.tol, options);
else
    [coords, e] = eigs(G.L, sparse(diag(G.d)), dim + 1, param.tol, options);
end
e = diag(e);

if nnz(e < param.tol) > 1
    warning('Multiple zero eigenvalues')
end

[~, ind] = sort(e, 'ascend');
coords = coords(:, ind);
coords = coords(:, 2:end);
end


