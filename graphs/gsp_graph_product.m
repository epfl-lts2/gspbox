function G = gsp_graph_product(G1,G2,param)
%GSP_GRAPH_PRODUCT Compute graph product between G1 and G2
%   Usage: G = gsp_graph_product(G1,G2, param)
%
%   Input parameters:
%       G1       : Graph 1
%       G2       : Graph 2
%       param    : Struct of parameters
%   Output parameters:
%       G        : Graph product
%
%   This function computes the graph product between two graphs according to rule specified in the parameters structure.
%   Factors graph structs are added to the struct of G.
%  
%   Additional parameters
%   ---------------------
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%    param.rule   : Graph product rule:
%     	+ cartesian (default)
%       + kronecker
%       + direct_sum
%       + strong
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_graph_product.php

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

% Author:  Francesco Grassi
% Date: July 2016

if nargin<3
    param = struct;
end

if ~isfield(param,'rule'), param.rule='cartesian';end

if ~isstruct(G1)
    W1 = G1;
else
    W1 = G1.W;   
end

if ~isstruct(G2)
    W2 = G2;
else
    W2 = G2.W;
end



switch param.rule
        
    case 'direct_sum'
        W = blkdiag(W1,W2);
        G = gsp_graph(W);
        
        C1=normc(G1.coords);
        C2=normc(G2.coords);
        C2(:,1)=C2(:,1)+max(abs(C1(:,1)))+.1;
        
        G.coords = [C1;C2];
        
        G.type = 'direct_sum';
        
    
    case 'cartesian'
        W = gsp_cartesian(W1,W2);
        G = gsp_graph(W);
        
        G.type = 'cartesian';

    case {'kronecker','direct'}
        W = kron(W1,W2);
        G = gsp_graph(W);
        G.type = 'kronecker';

    case 'strong'
        W = gsp_strong(W1,W2);
        G = gsp_graph(W);
        G.type = 'strong';   
        
    otherwise
        error('Unknown product rule. You can choose between Direct_Sum, Cartesian, Kronecker, Strong');
end


G.Gf = {G1,G2};


end


