function [ index ] = gsp_good_graph_index( G, X, param )
%GSP_GOOD_GRAPH_INDEX Index testing how well a given graph G, matches some data X
%   Usage: gsp_good_graph(G, X);
% 
%   Input parameters:
%       G          : the graph
%       X          : a data matrix
%       param      : structure of optional parameters
%   Output parameters: 
%       index      : the computed index 
% 
%   A wrapper function with which one may test how well a given graph G,
%   matches some data X.
% 
%   Example:::
%       G = gsp_2dgrid(16);
%       X = pinv(full(G.L)) * randn(G.N, G.N);
%       param.verbose = 1;
%       param.index = 'tcer';
%       index = gsp_good_graph_index(G, X, param)
%       param.index = 'stationarity';
%       index = gsp_good_graph_index(G, X, param)
%
%   Optional paramaters
%   -------------------
%
%   * *param.index*: 'tcer' or 'stationarity' (default 'tcer'). 
% 
%   See also: gsp_learn_tcer, gsp_stationarity_ratio


% Author  : Andreas Loukas
% Date    : 15 Nov 2016

% Handle input
if nargin < 3, param = struct(); end
if not(isfield(param, 'index'));   param.index = 'tcer'; end;

switch param.index,
  
    case 'tcer',    
        index = gsp_learn_tcer(G, X, param);
        
    case 'stationarity'        
        index = gsp_stationarity_ratio(G, X*X', param);
        
    otherwise, 
        error('uknown index.');

end
end

