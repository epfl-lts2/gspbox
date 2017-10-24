% sgwt_setpath : Set paths for SGWT toolbox
%
% user should set SGWT_ROOT to be directory where sgwt_toolbox is installed

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

SGWT_ROOT=[fileparts(mfilename('fullpath')),filesep];

sgwt_relpathlist={'',...
                  'utils',...
                  'demo',...
                  'mex',...
                 };
fprintf('Welcome to sgwt_toolbox. SGWT root directory is %s\n',SGWT_ROOT);
for k=1:numel(sgwt_relpathlist)
  sgwt_tmp_pathname=[SGWT_ROOT,sgwt_relpathlist{k}];
  fprintf('adding path %s\n',sgwt_tmp_pathname);
  addpath(sgwt_tmp_pathname);

end
clear sgwt_relpathlist sgwt_tmp_pathname
