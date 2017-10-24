function [ f ] = WGFT_igft(f_hat,V)
% Compute the inverse graph Fourier transform of a function f_hat on a 
% graph whose eigenvectors are the columns of V

    f=V*f_hat;

end

