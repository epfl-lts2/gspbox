function d = distanz(x,y,algtype)
% DISTANZ : calculates the distances between all vectors in x and y.
%
% usage:
%    d = distanz(x,y);
%
% inputs:
%    x      matrix with col vectors
%    y      matrix with col vectors (default == x)
%    algtype   the type of algorithm that is used (default==3)
%
% outputs:
%    d      distance matrix, not squared
%
% note:
%    part of the code is inspired by dist.m of the nntoolbox, other
%    part adapted from Francis Bach who took it from Roland
%    Bunschoten.
%
% sth * 19apr2002
% Adapted from create.m, originally written by
% (c) Stefan Harmeling, 2006

if exist('algtype')~=1 || isempty(algtype), algtype = 3; end
switch algtype
 case 1  % inspired by dist.m
  if exist('y')~=1 || isempty(y)
    % here comes code just for x
    [rx,cx] = size(x);
    d = zeros(cx,cx);
    nuller = zeros(cx,1);
    for c = 1:cx
      d(c,:) = sum((x-x(:,c+nuller)).^2,1);
    end
  else
    % here comes code for x and y
    [rx,cx] = size(x);
    [ry,cy] = size(y);
    if rx~=ry, error('x and y do not fit'), end
    d = zeros(cx,cy);
    if cx>cy
      nuller = zeros(cx,1);
      for c = 1:cy
	d(:,c) = sum((x-y(:,c+nuller)).^2,1)';
      end
    else
      nuller = zeros(cy,1);
      for c = 1:cx
	d(c,:) = sum((x(:,c+nuller)-y).^2,1);
      end
    end
  end
 
 case 2  % same as case 1, but with repmat instead of nuller
  if exist('y')~=1 || isempty(y)
    % here comes code just for x
    [rx,cx] = size(x);
    d = zeros(cx,cx);
    nuller = zeros(cx,1);
    for c = 1:cx
      d(c,:) = sum((x-repmat(x(:,c),[1 cx])).^2,1);
    end
  else
    % here comes code for x and y
    [rx,cx] = size(x);
    [ry,cy] = size(y);
    if rx~=ry, error('x and y do not fit'), end
    d = zeros(cx,cy);
    if cx>cy
      nuller = zeros(cx,1);
      for c = 1:cy
	d(:,c) = sum((x-repmat(y(:,c),[1 cx])).^2,1)';
      end
    else
      nuller = zeros(cy,1);
      for c = 1:cx
	d(c,:) = sum((repmat(x(:,c),[1 cy])-y).^2,1);
      end
    end
  end
  
 case 3  % inspired by Roland Bunschoten
  if exist('y')~=1 || isempty(y)
    % here comes code just for x
    cx = size(x,2);
    xx = sum(x.*x,1); xz = x'*x;
    d = abs(repmat(xx',[1 cx]) - 2*xz + repmat(xx,[cx 1]));
  else
    % here comes code for x and y
    [rx,cx] = size(x);
    [ry,cy] = size(y);
    if rx~=ry, error('x and y do not fit'), end
    xx = sum(x.*x,1); yy = sum(y.*y,1);  xy = x'*y;  
    d = abs(repmat(xx',[1 cy]) + repmat(yy,[cx 1]) - 2*xy);
  end
end

d = sqrt(d);

