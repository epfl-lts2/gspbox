% See http://stackoverflow.com/questions/16146599/create-artificial-data-in-matlab
% and http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
%
%   Url: http://lts2research.epfl.ch/gsp/doc/pointsclouds/to_be_included/datasetsdemo.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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

figure;
hold on;
dotsize = 12;
 colormap([1 0 .5;   % magenta
           0 0 .8;   % blue
           0 .6 0;   % dark green
           .3 1 0]); % bright green

subplot(231);
data = twospirals();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Two spirals');

subplot(232);
data = clusterincluster();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Cluster in cluster');

subplot(233);
data = corners();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Corners');

subplot(234);
data = halfkernel();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Half-kernel');

subplot(235);
data = crescentfullmoon();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Crescent & Full Moon');

subplot(236);
data = outlier();
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal;
title('Outlier');
