function [x,y] = gsp_twospirals(N, degrees, start, noise)
%GSP_TWOSPIRALS Generate "two spirals" dataset with N instances
%   Usage:  gsp_twospirals(N, degrees, start, noise);
%           gsp_twospirals(N, degrees, start);
%           gsp_twospirals(N, degrees);
%           gsp_twospirals(N);
%           gsp_twospirals();
%
%   Input arguments:
%           N       : Number of points (default 2000)
%           degrees : The length of the spirals in degree (default 570)
%           start   : how far from the origin the spirals start, in degrees (default 90)
%           noise   : Noise level (default 0.2) 
%
%   Output arguments:
%           x       : Position of the points
%           y       : label
%           
%   Note that for *noise=0*, there is no noise and at *noise=1* the spirals
%   will start overlapping 
%
%   This function is adaptated from:
%   http://stackoverflow.com/questions/16146599/create-artificial-data-in-matlab
%   and 
%   http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
%
%   Example:::
%
%       [x,y] = gsp_twospirals();
%       figure(); 
%       plot(x(y>0,1),x(y>0,2),'xr');
%       hold on
%       plot(x(y<=0,1),x(y<=0,2),'ob');
%       hold off
%       
%   See also: 

% Date: 15 August 2015
% Adaptated by: Nathanael Perraudin

    if nargin < 1
        N = 2000;
    end
    if nargin < 2
        degrees = 570;
    end
    if nargin < 3
        start = 90;
    end
    if nargin < 4
        noise = 0.2;
    end  
    
    deg2rad = (2*pi)/360;
    start = start * deg2rad;

    N1 = floor(N/2);
    N2 = N-N1;
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;   
    d1 = [-cos(n).*n + rand(N1,1)*noise sin(n).*n+rand(N1,1)*noise zeros(N1,1)];
    
    n = start + sqrt(rand(N1,1)) * degrees * deg2rad;      
    d2 = [cos(n).*n+rand(N2,1)*noise -sin(n).*n+rand(N2,1)*noise ones(N2,1)];
    
    data = [d1;d2];
    x = data(:,1:2);
    y = data(:,3);
end