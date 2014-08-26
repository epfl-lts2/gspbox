function [ G ] = gsp_create_laplacian( G,type )
%GSP_CREATE_LAPLACIAN create the graph laplacian of the graph G
%   Usage: G = gsp_create_laplacian( G,type );
%          G = gsp_create_laplacian( G );
%
%   Input parameters:
%       G   : Graph structure (or cell array of graph structure) 
%       type: Type of laplacian (string)
%   Output parameters:
%       G   : Graph structure (or cell array of graph structure) 
%
%   This function create the graph laplacian of the graph G and store it
%   into G.
%
%   The variable type contains the different laplacian type.
%
%    combinatorial*: Non normalized laplacian. This is the default.
%    normalized*: Normalized laplacian
%    none*: No laplacian
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_create_laplacian.php

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
% Date  : 09.12.2013

if numel(G)>1
    Ng = numel(G);
    for ii = 1:Ng
        if nargin<2
            G{ii} = gsp_create_laplacian(G{ii});
        else
            G{ii} = gsp_create_laplacian(G{ii}, type);
        end
    end     
    return;
end


if nargin<2
    if ~isfield(G,'lap_type')
        type='combinatorial';
        G.lap_type = type;
    else
        type = G.lap_type;
    end
end




switch type
    case 'combinatorial'
        G.L=diag(sum(G.W))-G.W;
    case 'normalized'
        D = diag(sum(G.W).^(-0.5));
        G.L=sparse(eye(G.N))-D*G.W*D;
    case 'none'
        G.L=sparse(0);
    otherwise
        error(' Unknown laplacian type')
end

G.lap_type = type;


end


