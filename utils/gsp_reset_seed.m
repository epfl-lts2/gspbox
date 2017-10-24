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

