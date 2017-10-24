function [G]=gsp_two_moons(sigmag,N,sigmad,d)
%GSP_TWO_MOONS  Initialize a two moons graph
%   Usage:  G = gsp_two_moons(sigmag);
%           G = gsp_two_moons(sigmag,N,sigmad);
%           G = gsp_two_moons(sigmag,N,sigmad,d);
%
%   Input parameters:
%         sigmag: parameter for the construction of the graph (default 0.05)
%         N     : Number of vertices. (default 2000)
%         sigmad: variance of the data 
%         d     : distance of the two moons (default 0.5)
%
%   Output parameters:
%         G     : Graph structure.
%
%   Initialise the two moons graph. If less than 2 parameters are given, a standard
%   graph is created. Otherwise, ....
%
%   Example:::
%
%          G = gsp_two_moons();
%          gsp_plot_signal(G,G.labels);
%

% Author :  Nathanael Perraudin

if nargin < 1
    sigmag = 0.05;
end

paramnn.sigma = sigmag;
paramnn.k = 5; 

if nargin < 2
    P = gsp_pointcloud('two_moons');
    G = gsp_nn_graph(P,paramnn);
    G.labels = (2*((1:G.N)>1000)-1)';
    G.type = 'Two moons standard';
else
        
    if nargin<4
       d = 0.5; 
    end
    N1 = floor(N/2);
    N2 = N - N1;

    % moon 1
    phi1 = rand(N1,1)  * pi;
    r1 = 1;
    rb = sigmad * randn(N1,1);
    %rb = laprnd(N1, 1, 0, sigmad);
    ab = rand(N1,1) * 2 * pi;
    b = rb .*exp(1i.*ab);
    bx = real(b);
    by = imag(b);
    moon1x = cos(phi1) .* r1 + 0.5 + bx;
    moon1y = -sin(phi1) .* r1 + by-(d-1)/2;
%     plot(moon1x , moon1y,'xb');

    % moon 2
    phi2 = rand(N2,1)  * pi;
    r2 = 1;
    rb = sigmad * randn(N2,1);
    ab = rand(N2,1) * 2 * pi;
    b = rb .*exp(1i.*ab);
    bx = real(b);
    by = imag(b);
    moon2x = cos(phi2) .* r2 - 0.5 + bx;
    moon2y = sin(phi2) .* r2 + by+(d-1)/2;
%     hold on;
%     plot(moon2x , moon2y,'xr');
    data = [moon1x,moon1y;moon2x,moon2y]';
    G = gsp_nn_graph(data',paramnn);
    G.labels = (2*((1:G.N)>N1)-1)';
    G.type = 'Two moons synthetised';
end



 




end
