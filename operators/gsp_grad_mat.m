function [ D ] = gsp_grad_mat( G )
%GSP_GRAD_MAT Gradient sparse matrix of the graph G
%   Usage:  D = gsp_gradient_mat(G);
%
%   Input parameters:
%       G   : Graph structure
%
%   Output parameters:
%       D   : Gradient sparse matrix
%
%   This function return the gradient matrix. To be more effiecient, call
%   the function: 
%
%           G = gsp_adj2vec(G)
%
%   before this function.
%
%   Example:
%
%       N = 40;
%       G = gsp_sensor(N);
%       G = gsp_adj2vec(G);
%       D = gsp_grad_mat(G);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_grad_mat.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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
% Date  : 14 Mai 2014
% Testing: test_operators


if ~isfield(G,'v_in')
    G = gsp_adj2vec(G);
    warning(['GSP_GRADIENT_MAT: To be more efficient you should run: ',...
        'G = gsp_adj2vec(G); before using this proximal operator.']);
end

if isfield(G,'Diff');
    D = G.Diff;
    return;
end
    

n = G.Ne;
Dr = [1:n 1:n];
Dc(1:n) = G.v_in;
Dc(n+1:2*n) = G.v_out;
Dv(1:n) = sqrt(G.weights);
Dv(n+1:2*n) = -sqrt(G.weights);

D = sparse(Dr,Dc,Dv,n,G.N);

end


