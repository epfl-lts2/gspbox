function [ Gout ] = gsp_separate_graph( G )
%GSP_SEPARATE_GRAPH Separate the graph G into disconnected components
%   Usage: [ Gout ] = gsp_separate_graph( G )
%
%   Input parameters:
%       G   : Original graph
%
%   Output parameters:
%       Gout: Cell array of graphs
%

% Author: Nathanael Perraudin
% Date : Mai 2014

% TODO: DO this operation better
G = gsp_compute_fourier_basis(G);

NG = sum(G.E < 1e-10);

error('Not implemented');

Gout = 0;


end

