% sgwt_view_design : display filter design in spectral domain
%
% function sgwt_view_design(g,t,arange)
%
% This function graphs the input scaling function and wavelet
% kernels, indicates the wavelet scales by legend, and also shows
% the sum of squares G and corresponding frame bounds for the transform.
%
% Inputs :
% g - cell array of function handles for scaling function and
%   wavelet kernels
% t - array of wavelet scales corresponding to wavelet kernels in g
% arange - approximation range

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

function sgwt_view_design(g,t,arange)
x=linspace(arange(1),arange(2),1e3);
clf
hold on

J=numel(g)-1;
co=get(gca,'ColorOrder');
co=repmat(co,[J,1]);
G=0*x;
for n=0:J
    plot(x,g{1+n}(x),'Color',co(1+n,:));
    G=G+g{1+n}(x).^2;
end
plot(x,G,'k');
[A,B]=sgwt_framebounds(g,arange(1),arange(2));

hline(A,'m:');
hline(B,'g:');
leglabels{1}='h';
for j=1:J
    if ~isempty(t)
        
        leglabels{1+j}=sprintf('t=%.2f',t(j));
    else
        leglabels{1+j}='';
    end
end

leglabels{J+2}='G';
leglabels{J+3}='A';
leglabels{J+4}='B';

%set(gca,'Ytick',0:3);

legend(leglabels)
title(['Scaling function kernel h(x), Wavelet kernels g(t_j x), Sum ' ...
       'of Squares G, and Frame Bounds']);

function hline(y,varargin)
  xl=xlim;
  plot(xl,y*[1 1],varargin{:});
