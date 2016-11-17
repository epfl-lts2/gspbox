function gsp_reset_seed(n)
%GSP_RESET_SEED Reset the seed of the random number generator
%   Usage:  gsp_reset_seed(n);
%
%   Input parameters:
%       n   : seed
%   Ouptut parameters:
%       none
%
%   This function resets the seed
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_reset_seed.php

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

% Authors: Nathanael Perraudin, Vassilis Kalofolias
% Date  : 21 May 2014
% 
% global GLOBAL_rand
% 
% if GLOBAL_rand
% s = rng;

if nargin<1
    n = 0;
end

if verLessThan('matlab', '7.12.0')  % release 2011a has "rng"
    rand('twister',n); %#ok<RAND>
else
    rng(n,'twister');
end

end


