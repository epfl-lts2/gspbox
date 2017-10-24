function s = gsp_jtv_filter_analysis_simple(G,g,X,param)

if nargin<4
    param = struct;
end



Xhat = gsp_jft_simple(G,X);

s = zeros(G.N,G.jtv.T,numel(g),size(X,3));

filt = gsp_jtv_filter_evaluate_simple(G,g,G.e,param);

% Here we need to use repmat to avoid both for loops but this function
% is written for testing purposes only
for ii = 1:numel(g)
    for jj = 1:size(X,3)
        s(:,:,ii,jj) = gsp_ijft(G,conj(filt(:,:,ii)).*Xhat);
    end
end

if sum(abs(imag(s(:)))) < 1e-11 *norm(s(:))
    s = real(s);
end

end