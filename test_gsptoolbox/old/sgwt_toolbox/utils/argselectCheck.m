% argselectCheck : Check if control parameters are valid 
%
%  function argselectCheck(control_params,varargin_in)
%
%  Inputs: 
%  control_params and varargin_in are both cell arrays
%  that are lists of pairs 'name1',value1,'name2',value2,...
%
%  This function checks that every name in varargin_in is one of the names
%  in control_params. It will generate an error if not, otherwise it will 
%  do nothing
%
%  This is used at beginning of function to simulate keyword argument
%  passing. Typical usage is
%
%  argselectAssign(control_params);
%  argselectCheck(control_params,varargin);
%  argselectAssign(varargin);
%
%  where control_params is a cell list of variable,value pairs containing
%  the default parameter values.
%
% See also argselectAssign
%
% Author : David K. Hammond, EPFL LTS2
% Date : December, 2007
% Project : common utilities

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

function argselectCheck(control_params,varargin_in)
[param_names{1:(length(control_params)/2)}]=control_params{1:2:end};
% check that passed in control parameters are valid
for j=1:2:length(varargin_in)
  if(isempty(strmatch(varargin_in{j},param_names)))
    error(['Invalid control parameter : ',varargin_in{j},...
           ' Valid options are',sprintf('\n'),...
          sprintf('%s \n',param_names{:})]);
  end
end
