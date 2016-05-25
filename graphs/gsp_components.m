function [Idx, Gz, bl] = gsp_components(G)
%GSP_COMPONENTS cuts non connected graph into several connected ones
%   Usage: Idx             = gsp_components(G);
%          [Idx, Gz]       = gsp_components(G);
%          [Idx, Gz, bl] = gsp_components(G);
%
%   Input parameters
%         G         : Graph
%
%   Output parameters
%         Idx       : Index vector
%         Gz        : Cell of sub-graphs
%         bl      : 0 or 1 (Connected or Disconnected graph)
%
%   This function cuts the non connected graph G into smaller
%   connected sub-graphs that are elements of the cell Gz. In order to
%   access the ith sub-graph one must use Gz{i}. The vector Idx contains 
%   elements from 1 to the number of subgraphs. The indices of the
%   elements of a specific subgraph correspond to the indices of the
%   elements of Idx that contain the same value. In order to obtain the
%   coordenets of the first subgraph one can use 
%   G.coords(Idx==1, :); 
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_components.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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
% Date    : 24/11/2015


if gsp_check_connectivity(G);
    Idx = ones(1,G.N);
    if nargout > 1 , Gz{1} = G; end
    if nargout > 2 , bl = 1;  end
    
    disp('The Graph is connected')
    
else
    disp('The Graph is not connected')
    disp('Component separation has begun')
    
    % if the structure G contains Fourier basis
    if isfield(G,'U') && isfield(G,'e') 
        tol = 1e-6;
        num_comp = nnz(G.e(1:end) < max(G.e)*tol);
        Idx = zeros(1, G.N);
        n=0;
        % while there are unlabeled elements in the vector Idx (elements
        % that are still 0) repeat. 
        while nnz(Idx) ~= G.N
            n=n+1;
            temp = bsxfun(@plus,G.U(:,1:num_comp),-G.U(find(Idx==0, 1)...
                ,1:num_comp));
            dist = sum(temp.^2,2);
            Idx(dist<=tol^2) = n;
        end
        
    else
        n = 0;
%        [D, ~] = gsp_weight2distance(G);
        Idx = zeros(1, G.N);
        % while there are unlabeled elements in the vector Idx (elements
        % that are still 0) repeat.
        while nnz(Idx) ~= G.N
            n = n + 1; % Label of subgraph elements
            d = dijkstra(G.W>0, find(Idx==0, 1)); % find unlabeled elements
            Idx(d~=inf) = n;
        end
    end
    
    % if the user also wants the subgraphs to be returned
    if nargout > 1
        comp = unique(Idx);  % Graph components
        Gz = cell(length(comp),1);
        for ii=1:length(comp)
            ind_comp = find(Idx == comp(ii));
            Gz{ii}.W = sparse(G.W(ind_comp, ind_comp));
            Gz{ii}.coords = G.coords(ind_comp, :);
            Gz{ii}.lap_type = G.lap_type;
            Gz{ii}.type = G.type;
            Gz{ii} = gsp_graph_default_parameters(Gz{ii});
        end
    end
    if nargout > 2
        bl = 0;
    end
end
end
