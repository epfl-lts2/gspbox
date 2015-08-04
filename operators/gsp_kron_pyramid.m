function [ Gs ] = gsp_kron_pyramid( G, Nlevels, param)
%GSP_KRON_PYRAMID compute a pyramid of graphs using the kron reduction
%   Usage: Gs = gsp_kron_pyramid( G, Nlevels);
%          Gs = gsp_kron_pyramid( G, Nlevels, param);
%
%   Input parameters:
%       G       : Graph structure
%       Nlevels : Number of level of decomposition
%       param   : Optional structure of parameters
%
%   Output parameters:
%       Gs      : Cell array of graphs
%
%   This function compute a pyramid of graph based on the Kron reduction.
%   The indices are taken as the positive entry of the highest eigenvector.
%
%   param is a structure of optional parameters containing the following
%   fields:
%
%    lambda*: Stability parameter. It add self loop to the graph to give
%     the alorithm some stability (default: 0.025).
%    sparsify*: Sparsify the graph after the Kron reduction (default: 1).
%    epsilon*: Sparsification parameter if the sparsification is used
%     (default:  min(2/sqrt(G.N), 0.1) ).
%    filters*: A cell array of filter that will be used for the analysis
%     and sytheis operator. If only one filter is given, it will be used
%     for all levels. You may change that later on. Default 
%
%            h(x) = 0.5 / ( 0.5 + x)
%   
%   Example:
%
%             N = 256;
%             G = gsp_sensor(N);
%             Nlevel = 5;
% 
%             Gs = gsp_kron_pyramid(G, Nlevel);
% 
%             figure;
%             for ii = 1:numel(Gs)
%                 subplot(2,3,ii)
%                 gsp_plot_graph(Gs{ii})
%                 title(['Reduction level: ', num2str(ii-1)]);
%             end
%
%   See also: gsp_pyramid_analysis gsp_pyramid_synthesis gsp_pyramid_cell2coeff
%
%   Demo: gsp_demo_pyramid
%
%   References:
%     D. I. Shuman, M. J. Faraji, and P. Vandergheynst. A framework for
%     multiscale transforms on graphs. arXiv preprint arXiv:1308.4942, 2013.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_kron_pyramid.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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
% Date  : 23 July 2014
% Testing: test_operators



if nargin < 3
    param = struct;
end
   
if ~isfield(param,'lambda'), param.lambda = 0.025; end;
if ~isfield(param,'sparsify'), param.sparsify = 1; end;
if ~isfield(param,'epsilon'), param.epsilon = min(10/sqrt(G.N), .1); end;

if ~isfield(param,'filters')
    filters=cell(Nlevels,1);
    for i=1:Nlevels
        filters{i}=@(x) .5./(.5+x);
    end
elseif length(param.filters)==1
    filters=cell(Nlevels,1);
    for i=1:Nlevels
        filters{i}=param.filters;
    end
elseif length(param.filters)==Nlevels
    filters=param.filters;
else
    error('param.filters should be a cell array of length 1 or num_levels');
end


Gs = cell(Nlevels + 1, 1);
Gs{1} = G;



for ii=1:Nlevels

    L_reg = full(Gs{ii}.L)+param.lambda*eye(Gs{ii}.N);
    [V,~] = eigs(L_reg, 1);
    % Select the bigger group
    V = V >= 0;
    if sum(V) >= Gs{ii}.N/2
        ind = find(V);
    else
        ind = find(1-V);
    end
    
    
    if param.sparsify
        Gtemp = gsp_kron_reduction(Gs{ii}, ind);
        Gs{ii+1} = gsp_graph_sparsify(Gtemp, ...
            max(param.epsilon, 2/sqrt(Gs{ii}.N)) );
    else
        Gs{ii+1} = gsp_kron_reduction(Gs{ii}, ind);
    end
    Gs{ii+1}.pyramid.ind = ind;
    Gs{ii+1}.pyramid.green_kernel = @(x) 1 ./ ( param.lambda + x );   
    Gs{ii+1}.pyramid.filter = filters{ii};
    Gs{ii+1}.pyramid.level = ii;
    Gs{ii+1}.pyramid.K_reg = gsp_kron_reduction(L_reg,ind);
end


end

 

