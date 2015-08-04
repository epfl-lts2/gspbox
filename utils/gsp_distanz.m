function d = gsp_distanz(x,y,P)
%GSP_DISTANZ calculates the distances between all vectors in x and y
%   Usage: d = gsp_distanz(x,y);
%
%   Input parameters:
%       x   : matrix with col vectors
%       y   : matrix with col vectors (default == x)
%       P   : distance matrix (default Identity)
%   Output parameters:
%       d   : distance matrix, not squared
%
%   This code compute the following
%
%      d = ( (x-y)^T P (x-y) )^(0.5)
%
%   for all vectors in x an y!
%   
%   This code is not optimized for memory, but for speed because it uses no
%   loops.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_distanz.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.4.0
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




% Handle Input parameters
if nargin<1, error('Not enought inputs parameters!'); end
if nargin<2,  y=x; end

% Get the size
[rx,cx] = size(x);
[ry,cy] = size(y);

% Verify the size
if rx~=ry, error('The sizes of x and y do not fit!'), end

if nargin < 3 
    xx = sum(x.*x,1); % ||x||^2
    yy = sum(y.*y,1); % ||y||^2 
    xy = x'*y;        % <y,x>
    % \|x-y\|^2 = ||x||^2 +||y||^2 - 2 <y,x> 
    d = abs(repmat(xx',[1 cy]) + repmat(yy,[cx 1]) - 2*xy);
else
    
    [rp,rp2] = size(P);
    if rx~=rp, error('The sizes of x and P do not fit!'), end
    if rp2~=rp, error('P must be square!'), end
    
    xx = sum(x .* (P* x),1 ); % x^T P x
    yy = sum(y .* (P* y),1 ); % y^T P y 
    xy = x'*(P*y);        % x^T P y
    yx = y'*(P*x);        % x^T P y

    d = abs(repmat(xx',[1 cy]) + repmat(yy,[cx 1]) - xy-yx);
end

if sum(d(:)<0)
    warning('gsp_distanz: P is not semipositive or x is not real!')
end
    
% Take the square root
d = sqrt(d);


