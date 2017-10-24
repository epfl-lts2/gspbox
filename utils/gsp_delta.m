function [ d, do ] = gsp_delta( G,at )
%GSP_DELTA 



if isstruct(G)
    N = G.N;
else
    N = G;
end


Nat = numel(at);
d = zeros(N,Nat);
for ii = 1:Nat;
    d(at(ii),ii) = 1;
end

if nargout>1
    if isfield(G,'Gm')
        do = zeros(G.Gm.N-G.N,Nat);
    else
        error('This function returns two arguments only for oose.')
    end
end

end

