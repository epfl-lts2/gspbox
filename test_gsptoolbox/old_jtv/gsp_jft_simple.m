function Xhat = gsp_jft_simple(G,X)

[T] = size(X,2);

Xhat = fft(gsp_gft(G,X),[],2)/sqrt(T);

end