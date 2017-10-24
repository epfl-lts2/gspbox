function [ modulated ] = WGFT_gmod(g,k,V)

% Modulate a signal g on a graph whose eigenvectors are the columns of V
% kernel g is input in the vertex domain
% k is in 1:N

N=size(V,1);
dc=V(:,1);
factor=1./dc;
modulated=factor.*V(:,k).*g;

% modulated=sqrt(N)*V(:,k).*g; OLD CODE - ADDED FACTOR TO ADJUST FOR
% NORMALIZED
end

