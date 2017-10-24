function [ translated ] = WGFT_gtrans(g,i,V)

% Note: input kernel in the vertex domain

%Translate a signal g on a graph whose eigenvectors are the columns of V
N=size(V,1);
translated=sqrt(N)*WGFT_gconv(g,sgwt_delta(N,i),V); % extra sqrt(N) factor added here instead of convolution

end

