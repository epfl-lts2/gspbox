function r = rescale_center(x)
%RESCALE_CENTER Rescaling the dataset

N=size(x,2);
x=x-repmat(mean(x,2),[1,N]);
c=max(abs(x(:)));
r=x/c;


end

