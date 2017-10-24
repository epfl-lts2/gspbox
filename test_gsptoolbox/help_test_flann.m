data = single(rand(2,1e3));
test = single(rand(2,1));

idx = flann_build_index(data,struct('algorithm','kdtree','trees',1));
res = flann_search(idx, test, 50, struct('checks',-1));    
% -1 is FLANN_CHECKS_UNLIMITED
flann_free_index(idx);

rmax = sqrt(max(sum( bsxfun(@minus, data(:,res), test).^2)))
% rmax is the distance of the worst neighbor found
linres = find( (sum( bsxfun(@minus, data, test).^2)) <= rmax^2 );
% linres are results from linear scan up to rmax
length(linres)
figure(1)
% plot data and neighbors to see problem
plot(data(1,:),data(2,:),'kx')
hold all;
axis equal;
plot(test(1),test(2),'bx',data(1,res),data(2,res),'ro');
rectangle('position',[test'-[rmax,rmax] [2*rmax 2*rmax]],'curvature',[1 1]);
hold off;

%%
figure(2)
param.use_flann = 1;
param.k = 50;
param.flann_checks = 256;
[res, ~, Dist] = gsp_nn_distanz(data,test,param);

Dist2 = sqrt(sum( bsxfun(@minus, data(:,res), test).^2))';
rmax = sqrt(max(sum( bsxfun(@minus, data(:,res), test).^2)))
% rmax is the distance of the worst neighbor found
linres = find( (sum( bsxfun(@minus, data, test).^2)) <= rmax^2 );
% linres are results from linear scan up to rmax
length(linres)

% plot data and neighbors to see problem
plot(data(1,:),data(2,:),'kx')
hold all;
axis equal;
plot(test(1),test(2),'bx',data(1,res),data(2,res),'ro');
rectangle('position',[test'-[rmax,rmax] [2*rmax 2*rmax]],'curvature',[1 1]);
hold off;