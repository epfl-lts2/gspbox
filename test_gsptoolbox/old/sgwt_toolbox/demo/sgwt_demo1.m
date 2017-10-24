% sgwt_demo1 : SGWT for swiss roll data set
%
% This demo builds the SGWT for the swiss roll synthetic data set. It
% computes a set of scales adapted to the computed upper bound on the
% spectrum of the graph Laplacian, and displays the scaling function and
% the scaled wavelet kernels, as well as the corresponding frame bounds. It
% then computes the wavelets centered on a single vertex, and displays
% them. This essentally reproduces figure 3 from 
% Hammond,Vangergheynst, Gribonval 2010.

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

% Create swiss roll point cloud
function sgwt_demo1 
close all
fprintf('Welcome to SGWT demo #1\n');

npoints=500;
fprintf('Creating Swiss Roll point cloud with %g points\n',npoints);
dataparams=struct('n',npoints,'dataset',-1','noise',0,'state',0);
r=create_synthetic_dataset(dataparams);
x=rescale_center(r.x);


fprintf('Computing edge weights and graph Laplacian\n');
% Compute Weighted graph adjacency matrix, and graph Laplacian
d=distanz(x);
s=.1;
A=exp(-d.^2/(2*s^2)); 
L=full(sgwt_laplacian(A));

fprintf('Measuring largest eigenvalue, lmax = ');

%% Design filters for transform
Nscales=4;
lmax=sgwt_rough_lmax(L);
fprintf('%g\n',lmax);

fprintf('Designing transform in spectral domain\n'); 

designtype='abspline3';
%designtype='mexican_hat';
%designtype='meyer';
%designtype='simple_tf';

[g,t]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);

arange=[0 lmax];
%% Display filter design in spectral domain
figure
sgwt_view_design(g,t,arange);
ylim([0 3])
set(gcf,'position',[0 780,600,300])
%% Chebyshev polynomial approximation
m=50; % order of polynomial approximation
fprintf('Computing Chebyshev polynomials of order %g for fast transform \n',m);
for k=1:numel(g)
  c{k}=sgwt_cheby_coeff(g{k},m,m+1,arange);
end

%% compute transform of delta at one vertex
jcenter=32; % vertex to center wavelets to be shown
fprintf('Computing forward transform of delta at vertex %g\n',jcenter);
N=size(L,1);
d=sgwt_delta(N,jcenter);
% forward transform, using chebyshev approximation
wpall=sgwt_cheby_op(d,L,c,arange);

fprintf('Displaying wavelets\n');
if (exist('OCTAVE_VERSION','builtin')~=0)
    % settings for octave
    msize=5;
    plotchar='s';
else    
    % settings for native matlab
    msize=100;
    plotchar='.';
end


cp=[-1.4,-16.9,3.4]; % camera position
%% Visualize result

% show original point
ws=300;
figure;
xp=0; yp=ws+100;
set(gcf,'position',[xp,yp,ws-10,ws+10]);
scatter3(x(1,:),x(2,:),x(3,:),msize,d,plotchar);
set(gcf,'Colormap',[.5 .5 .5;1 0 0]);
clean_axes(cp);
title(sprintf('Vertex %g',jcenter));

% show wavelets
for n=1:Nscales+1
    wp=wpall{n};
    figure
    xp=mod(n,3)*(ws+10);
    yp=(1-floor((n)/3))*(ws+100);
    set(gcf,'position',[xp,yp,ws-10,ws+10]);
    scatter3(x(1,:),x(2,:),x(3,:),msize,wp,plotchar);
    colormap jet
    caxis([-1 1]*max(abs(wp)));
    clean_axes(cp);

    hcb=colorbar('location','north');
    cxt=get(hcb,'Xtick');
    cxt=[cxt(1),0,cxt(end)];
    set(hcb,'Xtick',cxt);
    cpos=get(hcb,'Position');
    cpos(4)=.02; % make colorbar thinner
    set(hcb,'Position',cpos);
    set(hcb,'Position',[.25 .91 .6 .02]);
    
    if n==1
      title('Scaling function');
    else      
      title(sprintf('Wavelet at scale j=%g, t_j = %0.2f',n-1,t((n-1))));

    end
end


function clean_axes(cp)
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
set(gca,'Xtick',[-1 0 1]);
set(gca,'Ytick',[-1 0 1]);
set(gca,'Ztick',[-1 0 1]);
axis square
set(gca,'CameraPosition',cp);

% rescale_center
% center input data at origin, then rescale so that all coordinates
% are between -1 and 1
% 
% x should be dxN
function r=rescale_center(x)
N=size(x,2);
d=size(x,1);
x=x-repmat(mean(x,2),[1,N]);
c=max(abs(x(:)));
r=x/c;
