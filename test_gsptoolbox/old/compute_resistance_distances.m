function [resistance_distances] = compute_resistance_distances(L)

% This code is slow and should be improved

pseudo=pinv(full(L));
resistance_distances=zeros(size(L));
N=size(L,1);
for i=1:N
    for j=1:N
        resistance_distances(i,j)=pseudo(i,i)+pseudo(j,j)-2*pseudo(i,j);
    end
end

