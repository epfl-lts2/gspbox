function [G, tt]=gsp_swiss_roll(N,rand_state,param)
%GSP_SWISS_ROLL Initialize a swiss roll graph
%   Usage:  G = gsp_swiss_roll(N,rand_state,param);
%
%   Input parameters:
%         N          : Number of vertices.
%         s          : sigma ( default: sqrt(2/N))
%         thresh     : threshold (default: 1e-6)
%         rand_state : rand seed (default: 0)
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_create_swiss_roll(N,s,thresh,rand_state)' initializes a graph
%   structure containing the swiss roll graph
%
%   Example:::
%
%          G = gsp_swiss_roll(200);
%          gsp_plot_graph(G);
%

% Author : David I Shuman, Nathanael Perraudin
% Test: test_graphs

if nargin<1
    N = 200;
end

if nargin<2
   rand_state = 0; 
end

gsp_reset_seed(rand_state);



a = 1;   % swiss roll goes from a*pi to b*pi
b = 4;   
y = rand(2,N);
% uniform distribution along the manifold (in data space)
tt = sqrt((b*b-a*a)*y(1,:)+a*a);
tt = pi*tt;
% now tt should go from a*pi to b*pi
height = y(2,:);
x = [tt.*cos(tt)/b^2; height; tt.*sin(tt)/b^2];

if nargin<3
    param = struct;
end

if ~isfield(param,'k'), param.k = 6; end
G = gsp_nn_graph(x',param);
G.map_coord = y';
G.type = 'Swiss Roll';

end


