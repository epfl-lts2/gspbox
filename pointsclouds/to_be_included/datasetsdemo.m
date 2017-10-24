% See http://stackoverflow.com/questions/16146599/create-artificial-data-in-matlab
% and http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly

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