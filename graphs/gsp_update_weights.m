function [ G ] = gsp_update_weights( G, W )
%GSP_UPDATE_WEIGHTS update weight matrix of the graph G
%   Usage: G = gsp_update_weight( G, W );
%
%   Input parameters
%       G   : Graph
%       W   : W new weight matrix
%   Output parameters
%       G   : Graph
% 
%   This function will update the weight of graph and recompute the
%   Laplacian the new Laplacian. It will also recompute the value *G.lmax*
%   if it is present in the graph.
%
%   For now this function does not recompute the Fourier basis!
%
%   Example::
%          
%          G.W = Wnew;
%          G = gsp_graph_default_parameters( G );
%
%          
%   See also: gsp_graph_default_parameters gsp_graph
%

% Author: Nathanael Perraudin
% Date  : 19.11.2014




G.W = W;

if isfield(G,'U')
    G = rmfield(G,'U');
end

if isfield(G,'e')
    G =rmfield(G,'e');
end

if isfield(G,'mu')
    G =rmfield(G,'mu');
end

G = gsp_graph_default_parameters(G);

if isfield(G,'lmax')
    G = rmfield(G,'lmax');
    G = gsp_estimate_lmax(G);
end


end

