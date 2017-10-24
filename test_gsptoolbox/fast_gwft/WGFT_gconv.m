function [ f_conv_g ] = WGFT_gconv(f,g,V)
% Compute the generalized convolution product of two signals f and g on a 
% graph whose eigenvectors are the columns of V

% To do: need to decide if sqrt(N) goes here or in translation, and if here, 
% generalize the normalization for normalized graph Laplacian
f_conv_g=WGFT_igft(WGFT_gft(f,V).*WGFT_gft(g,V),V);

end

