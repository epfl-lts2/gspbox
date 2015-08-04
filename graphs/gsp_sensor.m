function [ G ] = gsp_sensor( N,param)
%GSP_SENSOR Create a random sensor graph
%   Usage: G = gsp_sensor( N );
%          G = gsp_sensor( );
%          G = gsp_sensor( N,param );
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
%   param is an optional structure with the following field
%
%    param.Nc : Minimum number of connection (default 2)
%    param.regular*: Flag to fix the number of connections to Nc (default 0)
%    param.verbose*: display parameter - 0 no log - 1 display the errors (default 1)
%    param.N_try*: Number of attempts to create the graph (default 50)
%    param.distribute*: To distribute the points more evenly (default 0)
%    param.connected*: To force the graph to be connected (default 1)
%
%
%   Example:
%
%          G = gsp_sensor(300);
%          paramplot.show_edges = 1;
%          gsp_plot_graph(G,paramplot);
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/graphs/gsp_sensor.php

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



% Date: 6 june 2013
% Author: Nathanael Perraudin


if nargin< 2
   param=struct;
end
if nargin<1
   N=64; 
end

if ~isfield(param, 'Nc'), param.Nc = 2; end
if ~isfield(param, 'regular'), param.regular = 0; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'N_try'), param.N_try = 50; end
if ~isfield(param, 'distribute'), param.distribute = 0; end
if ~isfield(param, 'connected'), param.connected = 1; end


if param.connected
    for n=1:param.N_try


        [W, XCoords, YCoords] = create_weight_matrix(N,param);


        if gsp_check_connectivity_undirected(W)
            break;
        elseif n==param.N_try
            fprintf(' Warning! Graph is not connected\n');
        end

    end
else
    [W, XCoords, YCoords] = create_weight_matrix(N,param);
end


% Return the values
G.W=sparse(W);
G.W = gsp_symetrize(G.W);

G.plotting.limits=[0,1,0,1];
G.N=N;
G.coords=[XCoords,YCoords];
if param.regular
    G.type='regular sensor';
else
    G.type='sensor';
end

G.directed=0;
        
G = gsp_graph_default_parameters( G );
end



function [W, XCoords, YCoords] = create_weight_matrix(N,param)


    XCoords = zeros(N,1);
    YCoords = zeros(N,1);
    if param.distribute
        mdim=ceil(sqrt(N));
        ind=1;
        for ii=0:mdim-1
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


    % Compute the distanz between all the points


    target_dist_cutoff = 2*N^(-0.5);
    T = .6; 
    s = sqrt(-target_dist_cutoff^2/(2*log(T)));
    d = gsp_distanz([XCoords,YCoords]'); 
    W = exp(-d.^2/(2*s^2)); 

    W=W-diag(diag(W));  
    
    if param.regular
        W = get_nc_connection(W,param);
    else
        W2 = get_nc_connection(W,param);
        W(W<T) = 0; % Thresholding to have sparse matrix 
        W(W2>0) = W2((W2>0));
    end
    
end


function W = get_nc_connection(W,param)
    N = size(W,1);
    Wtemp = W;
    W = zeros(size(W)); % Start with an empty matrix
    for ii=1:N
        l=Wtemp(ii,:);
       for jj=1:param.Nc
          [val,ind] = max(l);
          W(ii,ind) = val;
          l(ind)=0;
       end
    end
    W=(W+W')/2;
        
end

