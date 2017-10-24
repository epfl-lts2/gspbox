function Gn = gsp_copy_graph_attributes(G,type,Gn)
%GSP_COPY_GRAPH_ATTRIBUTES copy the parameter of the graph
%   Usage: Gn = gsp_copy_graph_attributes(G);
%
%   Input arguments:
%       G       : Graph structure
%       type    : flag to select what to copy (default 1)
%       Gn      : Graph structure (optional)
%
%   Output arguments:
%       Gn      : Partial graph structure
%
%   This function copy optional argument of a graph but not the graph
%   itself. If a graph is given as a third argument, then it will copy
%   everything in this graph. Otherwise, the function will create a new
%   empty graph.
%
%
%   The flag type can take the following value:
%
%   * 0: copy only the the graph structure
%       * lap_type
%       * plotting (but not plotting.limits)
%
%   * 1: copy additionaly the field
%       * plotting.limits
%       * coords
%



% Author: Nathanael Perraudin
% Date  : 24 July 2014


if nargin<2
    type = 1;
end

if nargin<3
    Gn = struct;
end

if isfield(G, 'lap_type'), Gn.lap_type = G.lap_type; end
if isfield(G, 'plotting'), Gn.plotting = G.plotting; end



if type
    if isfield(G, 'coords'), Gn.coords = G.coords; end
else
    if isfield(Gn.plotting,'limits')
        Gn.plotting = rmfield(Gn.plotting,'limits');
    end
end


end
