function [ f_hat ] = WGFT_gft(f,V)
% Compute the graph Fourier transform of a function f on a graph whose
% eigenvectors are the columns of V

f_hat=V'*f;

end

