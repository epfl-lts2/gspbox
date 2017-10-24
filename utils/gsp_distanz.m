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
%   .. D = ( (X-Y)^T P (X-Y) )^(0.5)
%
%   .. math:: D = \sqrt{(X-Y)^T P (X-Y)}
%
%   for all vectors in X an Y!
%   
%   This code is not optimized for memory, but for speed because it uses no
%   loops.
%

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


