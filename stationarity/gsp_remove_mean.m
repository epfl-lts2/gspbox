function [ X, mX ] = gsp_remove_mean( X, dim )
%GSP_REMOVE_MEAN Remove the mean of the signal X
%   Usage:  X  = gsp_remove_mean( X );
%           X  = gsp_remove_mean( X, dim );
%           Ç X, mX ]  = gsp_remove_mean( ... );
%
%   Input parameters:
%       X       : data
%       dim     : dimension (default 1)
%   Output parameters:
%       X       : centered data
%       mX      : mean
%   
%   This function remove efficiently the mean of X across line or column of
%   the data.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/stationarity/gsp_remove_mean.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

% Author: Nathanael Perraudin
% Date  : 24 July 2016
% Testing: test_gsp_remove_mean



if nargin<2
    dim = 1;
end

mX = mean(X,dim);

X = bsxfun(@minus,X,mX); 


end


