% sgwt_laplacian :  Compute graph laplacian from connectivity matrix
%
% function L = sgwt_laplacian(A,varargin)
%
% Connectivity matrix A must be symmetric.  A may have arbitrary
% non-negative values, in which case the graph is a weighted
% graph. The weighted graph laplacian follows the definition in
% "Spectral Graph Theory" by Fan R. K. Chung
%
% Inputs :
% A - adjacency matrix
% Selectable Input parameters :
% 'opt' - may be 'raw' or 'normalized' (default raw) to select
%         un-normalized or normalized laplacian
%
% Outputs : 
% L - graph Laplacian

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

function L = sgwt_laplacian(A,varargin)
  control_params={'opt','raw'}; % or normalized
  argselectAssign(control_params);
  argselectCheck(control_params,varargin);
  argselectAssign(varargin);
  
  N=size(A,1);
  if N~=size(A,2)
    error('A must be square');
  end
  degrees=vec(full(sum(A)));
  % to deal with loops, must extract diagonal part of A
  diagw=diag(A);
  
  % w will consist of non-diagonal entries only
  [ni2,nj2,w2]=find(A);
  ndind=find(ni2~=nj2); % as assured here
  ni=ni2(ndind);
  nj=nj2(ndind);
  w=w2(ndind);
  
  di=vec(1:N); % diagonal indices  
  
  switch opt
   case 'raw'
    % non-normalized laplacian L=D-A
    L=sparse([ni;di],[nj;di],[-w;degrees-diagw],N,N);
   case 'normalized'
    % normalized laplacian D^(-1/2)*(D-A)*D^(-1/2)
    % diagonal entries
    dL=(1-diagw./degrees); % will produce NaN for degrees==0 locations
    dL(degrees==0)=0;% which will be fixed here
    % nondiagonal entries
    ndL=-w./vec( sqrt(degrees(ni).*degrees(nj)) );
    L=sparse([ni;di],[nj;di],[ndL;dL],N,N);
   otherwise
    error('unknown option');
  end
