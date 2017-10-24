function bool = gsp_check_jtv(G)
%GSP_CHECK_JTV Check if G is a JTV graph
%   Usage:  bool = gsp_check_jtv(G):
%
%   Input parameters:
%       G           : Graph structure
%   Output parameters:
%       bool        : boolean
%
%   This function check if the structure G is a valid Joint Time-Vertex
%   Graph structure
%

bool = 0;

if ~isstruct(G)
    error('Input is not a valid Graph structure')
end

if ~isfield(G.jtv,'T')
    return
end

if or(~isfield(G.jtv,'T'),~isfield(G.jtv,'fs'))
    return
end

bool = 1;