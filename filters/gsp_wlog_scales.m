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

% Author: David K. Hammond, Nathanael Perraudin
% Date  : 18 March 2014


    t1=1;
    t2=2;

    smin=t1/lmax;
    smax=t2/lmin;
    % scales should be decreasing ... higher j should give larger s
    s=exp(linspace(log(smax),log(smin),Nscales));
  
end
