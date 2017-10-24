% test_cheby_op_mex : compare compiled (mex) and raw matlab versions of SGWT forward
% and adoint transformations

N=32^2;
A=sgwt_meshmat([sqrt(N),sqrt(N)]);

L=sgwt_laplacian(A);
lmax=sgwt_rough_lmax(L);
lpfactor=200;
nscales=5;
g=sgwt_filter_design(lmax,nscales,'lpfactor',lpfactor);

arange=[0,lmax];
m=50;
for j=1:nscales+1
    c{j}=sgwt_cheby_coeff(g{j},m,m+1,arange);
end

d=zeros(N,1);
d(ceil(N/2))=1;

% test forward transform
tic;y1=sgwt_cheby_op(d,L,c,arange);t1=toc;
tic;y2=cheby_op_mex(d,L,c,arange);t2=toc;
tic;y3=cheby_op_mex(d);t3=toc;

for ks=1:nscales+1;
    forward_equal(ks)=isequal(y1{ks},y2((1:N)+(ks-1)*N),y3((1:N)+(ks-1)*N) );
end

fprintf('all scales of forard transform equal? %d \n',all(forward_equal));
% test adjoint
tic;adj1=sgwt_adjoint(y1,L,c,arange);t1adj=toc;
tic;adj2=cheby_op_adjoint_mex(y2,L,c,arange);t2adj=toc;
tic;adj3=cheby_op_adjoint_mex(y2);t3adj=toc;

adjequal=isequal(adj1,adj2,adj3);
fprintf('adjoint transform equal? %d \n',adjequal);
fprintf('forward runtimes : raw matlab %g, mex %g (speedup %g)\n',t1,t2,t1/t2);
fprintf('adjoint runtimes : raw matlab %g, mex %g (speedup %g)\n',t1adj,t2adj,t1adj/t2adj);

fprintf('forward runtimes : raw matlab %g, primed mex %g (speedup %g)\n',t1,t3,t1/t3);
fprintf('adjoint runtimes : raw matlab %g, primed mex %g (speedup %g)\n',t1adj,t3adj,t1adj/t3adj);
