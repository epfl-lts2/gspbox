% sgwt_show_im : Display image, with correct pixel zoom 
%
% sgwt_show_im(im,range,zoom)
%
% Inputs :
% im - 2-d image
% range - 2 element vector giving display color map range, 
% range(1) maps to black, range(2) maps to white
% If range not given, or empty matrix given for range, then
% the default is to set it to the minimum and maximum of input image.
% zoom - # of screen pixels taken by single image pixel. Default is 1

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

function sgwt_show_im(im,range,zoom)
  if nargin<3
    zoom=1;
  end
  if ( nargin<2 || isempty(range) )
    range(1)=min(im(:));
    range(2)=max(im(:));
  end
  
  nshades=256;
  d_im = ( im-range(1) ) *(nshades-1) /(range(2)-range(1));
  dsize=size(im)*zoom; % size in pixels to show on screen
  image( d_im );   
  colormap(gray(nshades));
  
  ax=gca;
  oldunits=get(ax,'Units');
  set(ax,'Units','pixels');
  pos = get(ax,'Position');
  axis('off');
  ctr = pos(1:2)+pos(3:4)/2;
  set(ax,'Position',[floor(ctr-dsize/2)+0.5, dsize] );
  axis('equal'); 

  % restore units
  set(ax,'Units',oldunits);
