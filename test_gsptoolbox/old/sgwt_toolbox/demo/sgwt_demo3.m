% sgwt_demo3 : Image decomposition with SGWT wavelets based on local adjacency.
%
% This demo builds the SGWT transform on a graph representing 
% adjacency on a pixel mesh with 4-nearest neighbor connectivity.
% This demonstrates inverse on problem with large dimension.
%
% The demo loads an image file and decomposes the image with the SGWT,
% showing the coefficients as images at each scale. The demo does not show
% the individual wavelets (this could be done by replacing the input 
% image by a "delta image" with a single unit nonzero pixel) .
%
% The inverse is then computed, from the original coefficients as well as 
% from a modified set of coefficients where only coefficients at one
% scale are preserved. This shows that the SGWT can generate a
% multiresolution decomposition for images. We don't claim that this
% particular local-adjacency based transform is better for image
% processing than other available wavelet image decompositions, but it
% demonstrates the flexibility of the SGWT.  

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

function sgwt_demo3
global SGWT_ROOT
close all;
fprintf('Welcome to SGWT demo #3\n');
% load image
imname = [fileparts(mfilename('fullpath')),filesep,'paques_attack.png'];

fprintf('loading image %s\n',imname);
im = double( imread(imname) );
% build mesh adjacency graph
fprintf('Building mesh adjacency graph\n');
A=sgwt_meshmat(size(im));
% transform
fprintf('Calculating graph Laplacian\n');
L=sgwt_laplacian(A);
fprintf('Measuring largest eigenvalue, lmax = ');
lmax=sgwt_rough_lmax(L);
arange=[0,lmax];
fprintf('%g\n',lmax);

Nscales=5;
fprintf('Designing transform in spectral domain\n');
[g,t]=sgwt_filter_design(lmax,Nscales);

m=25; % order of polynomial approximation
fprintf('Computing Chebyshev polynomials of order %g for fast transform \n',m);
for k=1:numel(g)
    c{k}=sgwt_cheby_coeff(g{k},m,m+1,arange);
end

fprintf('Computing forward transform\n');
wpall=sgwt_cheby_op(im(:),L,c,arange);

% invert with all subbands
fprintf('Computing inverse transform with all coefficients\n');
imr1=sgwt_inverse(wpall,L,c,arange);
imr1=reshape(imr1,size(im));

ks=3; % scale at which to keep coefficients, set all others to zero.
fprintf('\nsetting all coefficients to zero except wavelet scale %g\n',ks-1);
% invert with only one scale
for k=1:numel(wpall)
    wpall2{k}=zeros(size(wpall{k}));
end
wpall2{ks}=wpall{ks};
fprintf('Computing inverse transform with coefficients from wavelet scale %g only\n',ks-1);
imr2=sgwt_inverse(wpall2,L,c,arange);
imr2=reshape(imr2,size(im));

%% display results
figure(1)
set(gcf,'position',[ 5   730   350   350]);
sgwt_show_im(im)
title('original image');
set(gcf,'menubar','none')
figure(2)
set(gcf,'position',[365 730 350 350]);
sgwt_show_im(imr1)
title('reconstuction from all coefficients');
set(gcf,'menubar','none')

figure(3)
set(gcf,'position',[725 730 350 350]);
sgwt_show_im(imr2);
title(sprintf('reconstruction only from wavelets at scale %g',ks-1));
set(gcf,'menubar','none')

figure(4)
set(gcf,'position',[0 0 1150 700]);
set(gcf,'menubar','none')
for k=1:Nscales+1
    subplot(2,3,k);
    sgwt_show_im(reshape(wpall{k},size(im)));
    if k==1
        title('Scaling function coefficients');
    else
        title(sprintf('Wavelet coefficients at scale %g',k-1));
    end
end
