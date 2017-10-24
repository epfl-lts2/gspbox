function G = gsp_sensor(N, param)
%GSP_SENSOR Create a random sensor graph
%   Usage: G = gsp_sensor(N);
%          G = gsp_sensor( );
%          G = gsp_sensor(N, param);
%
%   Input parameters
%       - N     : Number of nodes (default 128)
%       - param : Structure of optional parameters
%   Output parameters
%       - G     : Graph
%
%   This function creates a 2 dimensional random sensor graph. All the
%   coordinates are between 0 and 1.
%
%   *param* is an optional structure with the following field
%
%   * *param.verbose*: display parameter - 0 no log - 1 display the errors (default 1)
%   * *param.N_try*: Number of attempts to create the graph (default 10)
%   * *param.distribute*: To distribute the points more evenly (default 0)
%   * *param.connected*: To force the graph to be connected (default 1)
%   * *param.nnparam*: optional parameter for the gsp_nn_graph
%
%
%   Example:::
%
%          G = gsp_sensor(300);
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%



% Date: 6 june 2013
% Author: Nathanael Perraudin


if nargin < 2
    param = struct;
end
if nargin < 1
    N = 64;
end

if N<6
    error('N needs to be greater than 6')
end

if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'N_try'), param.N_try = 10; end
if ~isfield(param, 'distribute'), param.distribute = 0; end
if ~isfield(param, 'connected'), param.connected = 1; end
if ~isfield(param, 'nnparam'), param.nnparam = {}; end
if ~isfield(param.nnparam, 'k'), param.nnparam.k = 6; end



if param.connected
    for n=1:param.N_try
        [ XCoords, YCoords] = create_coords(N,param.distribute);
        % sort rows for plotting reasons
        G = gsp_nn_graph(sortrows([XCoords, YCoords]),param.nnparam);
        if gsp_check_connectivity(G)
            break;
        elseif n == param.N_try
            fprintf('Warning! Graph is not connected\n');
        end
    end
else
    [ XCoords, YCoords] = create_coords(N,param.distribute);
    % sort rows for plotting reasons
    G = gsp_nn_graph(sortrows([XCoords, YCoords]),param.nnparam);
end

% Return the values

G.type = 'sensor';


G = gsp_graph_default_parameters(G);
end



function [ XCoords, YCoords] = create_coords(N, distribute)

% TODO: VECTORIZE!!!!
XCoords = zeros(N,1);
YCoords = zeros(N,1);
if distribute
    mdim = ceil(sqrt(N));
    ind = 1;
    for ii = 0:mdim-1
        for jj=0:mdim-1
            if ind<=N
                XCoords(ind) = 1/mdim*rand(1)+ii*1/mdim;
                YCoords(ind) = 1/mdim*rand(1)+jj*1/mdim;
            end
            ind = ind+1;
        end
    end
else
    % take random coordinates in a 1 by 1 square
    XCoords = rand(N,1);
    YCoords = rand(N,1);
end


end

