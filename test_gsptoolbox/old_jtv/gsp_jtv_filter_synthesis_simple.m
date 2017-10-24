function X = gsp_jtv_filter_synthesis_simple(G,g,S,param)

if nargin<4
    param = struct;
end



X = zeros(G.N,G.jtv.T,size(S,4));
filt = gsp_jtv_filter_evaluate_simple(G,g,G.e,param);

% using repmat we should be able to get rid of one for loop
for jj = 1:size(S,4)
    for ii = 1:length(g)
        Shat = gsp_jft(G,S(:,:,ii,jj));
        X(:,:,jj) = X(:,:,jj) + gsp_ijft(G,filt(:,:,ii).*Shat);
    end
end

if sum(abs(imag(X(:)))) < 1e-11 *norm(X(:))
    X = real(X);
end

end