% sgwt_meshmat : Adjacency matrix for regular 2d mesh 
%
% function A=meshmat_p(dim,varargin)
%
% Inputs:
% dim - size of 2d mesh 
% Selectable control parameters:
% boundary - 'rectangle' or 'torus'
%
% Outputs:
% A - adjacency matrix

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

function A=sgwt_meshmat(dim,varargin)
  control_params={'boundary','rectangle'};
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);
  if (numel(dim)==1)
    dim=[1 1]*dim;
  end
  % build adjacency matrix : find i,j coordinates of center points
  % and right and bottom neighbors, then build connectivity matrix.
  % For each valid center,neighbor pair, will add A(center,neighbor)=1
  % and A(neighbor,center)=1, so A will be symmetric
  N=prod(dim);
  [alli,allj]=find(ones(dim));
  % (ci(k),cj(k)) has neighbor (ni(k),nj(k))
  ci=[alli;alli];
  cj=[allj;allj];
  ni=[alli  ; alli+1];
  nj=[allj+1; allj];
  switch boundary
   case 'rectangle' 
    % prune edges at boundary
    valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2));
    ni=ni(valid);
    nj=nj(valid);
    ci=ci(valid);
    cj=cj(valid);
    cind=dim(1)*(cj-1)+ci;
    nind=dim(1)*(nj-1)+ni;
   case 'torus'
    % wrap indices to make torus
    ni=mod(ni,dim(1))+1;
    nj=mod(nj,dim(2))+1;    
    cind=dim(1)*(cj-1)+ci;
    nind=dim(1)*(nj-1)+ni;
   otherwise
    error('unknown boundary option');
  end
  % assemble connection matrix
  A=sparse([cind,nind],[nind,cind],ones(1,2*numel(ni)),N,N);
  
