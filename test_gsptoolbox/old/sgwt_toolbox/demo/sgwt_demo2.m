% sgwt_demo2 : Allows exploring wavelet scale and approximation accuracy
%
% This demo builds the SGWT for the minnesota traffic graph, a graph
% representing the connectivity of the minnesota highway system. One center
% vertex is chosen, and then the exact (naive forward transform) and the
% approximate (via chebyshev polynomial approximation) wavelet transforms
% are computed for a particular value of the wavelet scale t. The relative
% error of the exact and approximate wavelets is computed. The user may
% then adjust the value of t, the degree m of the chebyshev polynomial
% approximation, and the center vertex in order to explore their effects.

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

function sgwt_demo2
close all
fprintf('Welcome to SGWT demo #2\n');

% touch variables to be shared among sub-functions
gb=[]; c=[]; 

% create UI elements
fh=figure('Visible','on','Name','demo 2 ui','Position',[425,920,400,150]);
uipanelh=uipanel('Parent',fh,'Title','','Units','pixels','BorderType','none');
tsliderh=uicontrol(uipanelh,'style','slider','max',50,'min',0,'value',1,...
                   'sliderstep',[.005 .1],'position',[25,10,300,20],...
                   'callback',{@tslider_callback});
msliderh=uicontrol(uipanelh,'style','slider','max',100,'min',1,'value',20,...
                   'sliderstep',[.001 .1],'position',[25,60,300,20],...
                   'callback',{@mslider_callback});
jbuttonh=uicontrol(uipanelh,'style','pushbutton','position',[50,110,150,20],...
                   'string','Select center vertex','callback',{@jbutton_callback});
ttexth=uicontrol(uipanelh,'style','text','string','','position',[325,10,100,20]);
mtexth=uicontrol(uipanelh,'style','text','string','','position',[325,60,100,20]);
jtexth=uicontrol(uipanelh,'style','text','string','','position',[325,100,100,20]);
uicontrol(uipanelh,'style','text','string',...
          'Chebyshev polynomial order (m)','position',...
           [60,80,200,20]);
uicontrol(uipanelh,'style','text','string',...
          'Wavelet scale (t)','position',...
          [60,30,200,20]);


%% Load graph and compute Laplacian
fprintf('Loading minnesota traffic graph\n');
Q=load('minnesota.mat');
xy=Q.xy;
A=Q.A;
N=size(A,1);
x=xy(:,1);
y=xy(:,2);

fprintf('Computing graph laplacian\n')
[ki,kj]=find(A);
L=sgwt_laplacian(A);
fprintf('Measuring largest eigenvalue, lmax = ');
lmax=sgwt_rough_lmax(L);
fprintf('%g\n',lmax);
arange=[0 lmax];

msize=100;

% initial values
t=3; % wavelet scale

m=20; % chebyshev polynomial order, for approximation
jcenter=550;

fprintf('\n');
update_uitext;
update_graphfig
update_kernel
update_waveletfigs

function update_graphfig
  figure(2)
  set(gcf,'renderer','zbuffer');
  fprintf('Displaying traffic graph\n');
  set(gcf,'position',[0,600,400,400]);
  %clf('reset');
  hold on
  scatter(x,y,msize,[.5 .5 .5],'.');
  plot([x(ki)';x(kj)'],[y(ki)';y(kj)'],'k');
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  axis equal
  axis off
  scatter(x(jcenter),y(jcenter),msize,'r.');
  drawnow
  end

  function update_kernel
  % select wavelet kernel
  t1=1;
  t2=2;
  a=2;
  b=2;
  tmin=t1/lmax; 
% scales t<tmin will show same wavelet shape as t=tmin, as 
% wavelet kernel g is monomial in interval [0,1)
set(tsliderh,'min',tmin);
  gb= @(x) sgwt_kernel_abspline3(x,a,b,t1,t2);
  g=@(x) gb(t*x);
  % polynomial approximation
  for k=1:numel(g)
    c=sgwt_cheby_coeff(g,m,m+1,arange);
  end
  lambda=linspace(0,lmax,1e3);
  figure(3)
  set(gcf,'position',[425,580,600,250])
  plot(lambda,g(lambda),lambda,sgwt_cheby_eval(lambda,c,arange));
  legend('Exact Wavelet kernel','Chebyshev polynomial approximation');
end

function update_waveletfigs
  
  fprintf('\nReomputing wavelets with t=%g, m=%g\n',t,m);
  d=sgwt_delta(N,jcenter);
  fprintf('Computing wavelet by naive forward transform\n');
  figure(4)
  set(gcf,'position',[0,100,400,400])
  wp_e=sgwt_ftsd(d,gb,t,L);
  show_wavelet(wp_e,x,y);
  % show wavelet (naive)
  title('exact wavelet (naive forward transform)');  
  fprintf('Computing wavelet by Chebyshev approximation\n');
  figure(5)
  set(gcf,'position',[425,100,400,400])
  % show wavelet (chebyshev)
  wp_c=sgwt_cheby_op(d,L,c,arange);
  show_wavelet(wp_c,x,y);
  title('approximate wavelet (transform via chebyshev approximation)');
  relerr=norm(wp_e-wp_c)/norm(wp_e);
  fprintf('Relative error between exact and approximate wavelet %g\n',relerr)
end

function show_wavelet(wp,x,y)
[Fs,s_ind]=sort(abs(wp),'descend');
scatter(x(s_ind),y(s_ind),msize,wp(s_ind),'.');
caxis([-1 1]*max(abs(wp)));
hcb=colorbar('location','north');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
cxt=get(hcb,'Xtick');
cxt=[cxt(1),0,cxt(end)];
set(hcb,'Xtick',cxt);
cpos=get(hcb,'Position');
cpos(4)=.02; % make colorbar thinner
set(hcb,'Position',cpos);
axis equal
axis off
end

function update_uitext
set(ttexth,'string',sprintf('t=%0.3f',t));
set(mtexth,'string',sprintf('m=%g',m));
set(jtexth,'string',sprintf('j=%g',jcenter));
end

function tslider_callback(source,eventdata)
t=get(tsliderh,'value');
update_uitext;
update_kernel;
update_waveletfigs;
end

function mslider_callback(source,eventdata)
newm=get(msliderh,'value');
if newm<m
  m=floor(newm);
else
  m=ceil(newm);
end
set(msliderh,'value',m);
update_uitext;
update_kernel;
update_waveletfigs;
end

function jbutton_callback(source,eventdata)
figure(2)
fprintf('Select new center vertex\n');
[xp,yp]=ginput(1);
oldjcenter=jcenter;
jcenter=argmin((xp-x).^2+(yp-y).^2);
scatter(x(jcenter),y(jcenter),msize,'r.');
scatter(x(oldjcenter),y(oldjcenter),msize,[.5 .5 .5],'.');
drawnow
update_uitext
update_waveletfigs
end

end

function i = argmin(x)
  i=min(find(x==min(x)));
end
