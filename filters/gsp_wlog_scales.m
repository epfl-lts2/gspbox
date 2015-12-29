function s = gsp_wlog_scales(lmin, lmax, Nscales)
%GSP_WLOG_SCALES compute logarithm scales for wavelet
%   Usage: s = gsp_wlog_scales(lmin, lmax, Nscales);
%
%   Input parameters:
%       lmin    : Minimum non zero eigenvalue
%       lmax    : Maximum eigenvalue
%       Nscales : Number of scale
%   Output parameters:
%       s       : scale
%
%   returns a (possibly good) set of wavelet scales given minimum nonzero
%   and maximum eigenvalues of laplacian
% 
%   returns scales logarithmicaly spaced between minimum and maximum
%   "effective" scales : i.e. scales below minumum or above maximum
%   will yield the same shape wavelet (due to homogoneity of kernel : 
%   currently assuming sgwt kernel g given as abspline with t1=1, t2=2)
%
%   Note that in design of transform with scaling function, lmin may be
%   taken just as a fixed fraction of lmax, and may not actually be the
%   smallest nonzero eigenvalue  
%   
%   This function is inspired by the sgwt_toolbox.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_wlog_scales.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Author: David K. Hammond, Nathanael Perraudin
% Date  : 18 March 2014


    t1=1;
    t2=2;

    smin=t1/lmax;
    smax=t2/lmin;
    % scales should be decreasing ... higher j should give larger s
    s=exp(linspace(log(smax),log(smin),Nscales));
  
end

