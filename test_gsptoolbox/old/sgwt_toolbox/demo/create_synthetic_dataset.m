function data = create_synthetic_dataset(data)
% create_synthetic_dataset creates test data for running nldr algorithms.
%
% inputs:
%    data      a struct describing the test data
%              .dataset the number of the example, see code for more infos
%              .n       the number of data points (default=400)
%              .state   the initial state for the random numbers (default=0)
%              .noise   the variance of Gaussian noise to add (default=0)
%              other options for some of the data sets (see code)
%              alternatively, data = 1 chooses the dataset directly,
%              the number of points defaults to 1000
%
% outputs:
%    data      a struct containing .x the generated data, each column is
%              a data point, and other stuff:
%              .z     the "correct" embedding
%              .e     some random noise of same dimensionality
%              .x_noisefree  the noisefree version of .x, i.e.
%                     .x = .xnoise_free + sqrt(.noise) * .e
%
% Adapted from create.m, originally written by
% (c) Stefan Harmeling, 2006
% using the examples of the original LLE and ISOMAP code.

if ~isfield(data, 'dataset'), 
  number = data;
  clear data
  data.dataset = number;
end
if ~isfield(data, 'n'), data.n = 400; end
if ~isfield(data, 'noise'), data.noise = 0.0; end
if ~isfield(data, 'state'), data.state = 0; end

% set the randomness
rand('state', data.state);
randn('state', data.state);

data.typ = 'data';
switch data.dataset
 case 0 % "swiss roll with hole" 
  data.name = 'swiss roll with hole';
  n = data.n;
  a = 1;   % swiss roll goes from a*pi to b*pi
  b = 4;   
  y = rand(2,n);
  % punch a rectangular hole at the center
  l1 = 0.05; l2 = 0.15;
  y = y - 0.5;
  ok = find((abs(y(1,:))>l1) | (abs(y(2,:))>l2));
  i = length(ok);
  y(:, 1:i) = y(:, ok);
  while (i<n)
    p = rand(2,1) - 0.5;
    if (abs(p(1))>l1) || (abs(p(2))>l2)
      i = i + 1;
      y(:,i) = p;
    end
  end
  y = y + 0.5;
  tt = (b-a)*y(1,:) + a;
  tt = pi*tt;
  height = 21*y(2,:);
  data.col = tt;
  data.x = [tt.*cos(tt); height; tt.*sin(tt)];
  data.z = [tt; height]; % the ground truth
  data.az = -4;
  data.el = 13;
  
 case -1 % "swiss roll" dataset extracted from LLE's swissroll.m
  data.name = 'uniform swiss roll';
  n = data.n;
  a = 1;   % swiss roll goes from a*pi to b*pi
  b = 4;   
  y = rand(2,n);
  data.z = y;  % the ground truth
  switch 1
   case 1
    % uniform distribution along the manifold (in data space)
    tt = sqrt((b*b-a*a)*y(1,:)+a*a);
   case 2
%    error('do not use this case')
    % nonuniform distribution along the manifold (in data space)
    tt = (b-a)*y(1,:) + a;  
  end
  tt = pi*tt;
  % now tt should go from a*pi to b*pi
  height = 21*y(2,:);
  data.col = tt;
  data.x = [tt.*cos(tt); height; tt.*sin(tt)];
  data.az = -4;
  data.el = 13;

 case 1 % "swiss roll (uniform in embedding space)" 
  % dataset extracted from LLE's swissroll.m
  data.name = 'classic swiss roll';
  n = data.n;
  a = 1;   % swiss roll goes from a*pi to b*pi
  b = 4;   
  y = rand(2,n);
  tt = (b-a)*y(1,:) + a;
  tt = pi*tt;
  height = 21*y(2,:);
  data.col = tt;
  data.x = [tt.*cos(tt); height; tt.*sin(tt)];
  data.z = [tt; height]; % the ground truth
  data.az = -4;
  data.el = 13;
  
 case 11 % "undersampled swiss roll"
  % dataset extracted from LLE's swissroll.m
  data.name = 'undersampled swiss roll';
  data.n = 100;
  n = data.n;
  a = 1;   % swiss roll goes from a*pi to b*pi
  b = 4;   
  y = rand(2,n);
  tt = (b-a)*y(1,:) + a;
  tt = pi*tt;
  height = 21*y(2,:);
  data.col = tt;
  data.x = [tt.*cos(tt); height; tt.*sin(tt)];
  data.z = [tt; height]; % the ground truth
  data.az = -4;
  data.el = 13;
  
 case 12 % "swiss roll"
  % dataset extracted from LLE's swissroll.m
  data.name = 'classic swiss roll';
  data.n = 400;
  n = data.n;
  a = 1;   % swiss roll goes from a*pi to b*pi
  b = 4;   
  y = rand(2,n);
  tt = (b-a)*y(1,:) + a;
  tt = pi*tt;
  height = 21*y(2,:);
  data.col = tt;
  data.x = [tt.*cos(tt); height; tt.*sin(tt)];
  data.z = [tt; height]; % the ground truth
  data.az = -4;
  data.el = 13;
  
 case 2 % "scurve" dataset extracted from LLE's scurve.m
  data.name = 'scurve';
  n = data.n;
  % I added 'ceil' and 'floor' to account for the case that n is odd
  angle = pi*(1.5*rand(1,ceil(n/2))-1); height = 5*rand(1,n);
  data.x = [[cos(angle), -cos(angle(1:floor(n/2)))]; height;[ sin(angle), 2-sin(angle)]];
  data.col = [angle, 1.5*pi + angle];
  data.z = [angle, 1.5*pi+angle; height]; % the ground truth
 
 case 3 % "square" dataset, a uniformly sampled 2D square randomly
         % rotated into higher dimensions
  data.name = 'square';
  n = data.n;
  d = 2;     % intrinsic dimension
  % optional parameter for dataset==3
  % data.D      dimension of the data
  if ~isfield(data, 'D'), data.D = 3; end
  % generate random rotation matrix
  D = data.D;
  A = randn(D, D);
  options.disp = 0;
  [R, dummy] = eigs(A*A', d, 'LM', options);
  tt = rand(d, n);
  data.col = tt(1,:);
  data.x = R*tt;
  data.z = tt;   % the ground truth
  data.az = 7;
  data.el = 40;
  
 case 4 % spiral: two dimensional "swiss roll"
  data.name = 'spiral';
  n = data.n;
  tt = (3*pi/2)*(1+2*rand(1, n));
  data.col = tt;
  data.x = [tt.*cos(tt); tt.*sin(tt)];
  data.z = tt; % the ground truth
  
 case -4 % spiral: two dimensional "swiss roll"
  data.name = 'noisy spiral';
  n = data.n;
  tt = (3*pi/2)*(1+2*rand(1, n));
  data.col = tt;
  data.x = [tt.*cos(tt); tt.*sin(tt)];
  data.x = data.x + randn(size(data.x));
  data.z = tt; % the ground truth
  
 case 5 % hole: a dataset with a hole
  data.name = 'hole';
  n = data.n;
  data.x = rand(2,n) - 0.5;
  % punch a rectangular hole at the center
  l1 = 0.2; l2 = 0.2;
  ok = find((abs(data.x(1,:))>l1) | (abs(data.x(2,:))>l2));
  i = length(ok);
  data.x(:, 1:i) = data.x(:, ok);
  while (i<n)
    p = rand(2,1) - 0.5;
    if (abs(p(1))>l1) || (abs(p(2))>l2)
      i = i + 1;
      data.x(:,i) = p;
    end
  end
  data.col = data.x(2,:);
  data.z = data.x;
  
 case 6 % P : taken from Saul's slides
  % note that for k=20, isomap and lle work fine which is very different
  % from the plots that Saul showed in his slides.
  data.name = 'P';
  load x
  x(2,:) = 500-x(2,:);
  data.x = x;
  data.z = x;
  data.col = data.z(2,:);
  data.n = size(x, 2);
  
 case 7 % fishbowl: uniform in data space
  gamma = 0.8;
  data.name = 'fishbowl (uniform in data space)';
  n = data.n;
  data.x = rand(3,n)-0.5;
  %project all data onto the surface of the unit sphere
  data.x = data.x ./ repmat(sqrt(sum(data.x.*data.x, 1)), [3 1]);
  ok = find(data.x(3,:) < gamma);
  i = length(ok);
  data.x(:, 1:i) = data.x(:, ok);
  while (i < n)
    p = rand(3,1)-0.5;
    p = p / sqrt(p'*p);
    if (p(3) < gamma)
      i = i+1;
      data.x(:, i) = p;
    end
  end
  % the projection on the plane works as follows:
  % start a beam from (0,0,1) through each surface point on the sphere
  % and look where it hits the xy plane.
  data.z = data.x(1:2,:) ./ repmat(1-data.x(3,:), [2 1]);
  data.col = data.x(3,:);
  data.az = -18;
  data.el = 16;
 case 8 % fishbowl: uniform in embedding space
  data.name = 'fishbowl (uniform in embedding space)';
  n = data.n;
  data.z = rand(2, n) - 0.5;
  % keep the disc
  ok = find(sum(data.z .* data.z) <= 0.25);
  i = length(ok);
  data.z(:, 1:i) = data.z(:, ok);
  while (i < n)
    p = rand(2,1) - 0.5;
    if (p'*p <= 0.25)
      i = i + 1;
      data.z(:, i) = p;
    end
  end
  gamma = 0.8;  % same role/parameter as in case 7
  data.z = 2*sqrt((1+gamma)/(1-gamma))*data.z;
  % project the disc onto the sphere
  alpha = 2 ./ (1 + sum(data.z .* data.z, 1));
  data.x = [repmat(alpha, [2 1]).*data.z; zeros(1, n)];
  data.x(3,:) = 1-alpha;
  data.col = data.x(3,:);
  data.az = -18;
  data.el = 16;
  
 case 9  % a gaussian blob
  data.name = 'gaussian blob';
  n = data.n;
  data.x = randn(3,n);
  data.z = data.x(2:3,:);
  data.col = data.x(3,:);
  
end


data.D = size(data.x, 1);  % dimensionality of the data
% finally generate noise
data.e = randn(size(data.x));
data.x_noisefree = data.x;  % the noise free data
data.x = data.x_noisefree + sqrt(data.noise)*data.e;

% precalculate the distanzmatrix
data.distances = distanz(data.x);


