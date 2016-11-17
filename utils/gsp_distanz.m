function D = gsp_distanz(X, Y, P)
%GSP_DISTANZ calculates the distances between all vectors in X and Y
%   Usage: D = gsp_distanz(X, Y);
%
%   Input parameters:
%       X   : matrix with col vectors
%       Y   : matrix with col vectors (default == X)
%       P   : distance matrix (default Identity)
%   Output parameters:
%       D   : distance matrix, not squared
%
%   This code computes the following
%
%      D = ( (X-Y)^T P (X-Y) )^(0.5)
%
%   for all vectors in X an Y!
%   
%   This code is not optimized for memory, but for speed because it uses no
%   loops.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_distanz.php

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

% Testing: test_gsp_distanz




% Handle Input parameters
if nargin<1, error('Not enought inputs parameters!'); end
if nargin<2,  Y=X; end

% Get the size
[rx, cx] = size(X);
[ry, cy] = size(Y);

% Verify the size
if rx~=ry, error('The sizes of x and y do not fit!'), end

if nargin < 3 
    xx = sum(X.*X,1); % ||x||^2
    yy = sum(Y.*Y,1); % ||y||^2 
    xy = X'*Y;        % <y,x>
    % \|x-y\|^2 = ||x||^2 +||y||^2 - 2 <y,x> 
    D = abs(repmat(xx',[1 cy]) + repmat(yy,[cx 1]) - 2*xy);
else
    
    [rp,rp2] = size(P);
    if rx~=rp, error('The sizes of x and P do not fit!'), end
    if rp2~=rp, error('P must be square!'), end
    
    xx = sum(X .* (P* X),1 ); % x^T P x
    yy = sum(Y .* (P* Y),1 ); % y^T P y 
    xy = X'*(P*Y);        % x^T P y
    yx = Y'*(P*X);        % x^T P y

    D = abs(repmat(xx',[1 cy]) + repmat(yy,[cx 1]) - xy-yx);
end

if sum(D(:)<0)
    warning('gsp_distanz: P is not semipositive or x is not real!')
end
    
% Take the square root
D = sqrt(D);

if nargin < 2   % if Y == X
    % The diagonal has to be zero!
    D(1:cx+1:end) = 0;
end



