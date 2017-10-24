function [f] = gsp_igft(G,f_hat)
%GSP_IGFT  Inverse graph Fourier transform
%   Usage:  f = gsp_igft(G,f_hat);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         f_hat      : Signal
%   Output parameters:
%         f          : Inverse graph Fourier transform of *f_hat*
%
%   'gsp_igft(G,f_hat)' computes a graph Fourier transform of the signal
%   *f_hat* with respect to the Fourier basis of the graph G: G.U.
%   Alternatively, one can provide directly the Fourier basis instead of
%   the graph G. 
%
%   .. f = U * f_hat 
%
%   .. math:: \hat{f}(\lambda_{\ell})=\langle f , u_{\ell} \rangle
%
%   To compute the Fourier basis of a graph G, you can use the function::
%
%           G = gsp_compute_fourier_basis(G);
%
%   Example:::
%
%       N = 30;
%       G = gsp_sensor(N);
%       G = gsp_compute_fourier_basis(G);
%       f_hat = zeros(N,1);
%       f_hat(5) = 1;
%       f = gsp_igft(G,f_hat);
%       gsp_plot_signal(G,f);  
%
%   See also: gsp_gft, gsp_gwft, gsp_compute_fourier_basis
%

% Author : Nathanael Perraudin, David I Shuma
% Date: 

s = size(f_hat);


if ~isnumeric(G);
    if ~gsp_check_fourier(G)
       error(['GSP_IGFT: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    % U = G.U;
    f=reshape(G.U * reshape(f_hat,G.N,[]),s);
else
    % U = G;
    f=reshape(G * reshape(f_hat,G.N,[]),s);
end


