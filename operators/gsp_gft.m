function [f_hat]=gsp_gft(G,f)
%GSP_GFT  Graph Fourier transform
%   Usage:  f_hat=gsp_gft(G,f);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         f          : f (signal)
%   Output parameters:
%         f_hat      : Graph Fourier transform of *f*
%
%   'gsp_gft(G,f)' computes a graph Fourier transform of the signal $f$
%   with respect to the Fourier basis of the graph G: G.U. Alternatively,
%   one can provide directly the Fourier basis instead of the graph G.
%
%   .. f_hat = U' * f 
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
%       f = sin((1:N)'*2*pi/N);
%       fhat = gsp_gft(G,f);
%       gsp_plot_signal_spectral(G,fhat);  
%
%   See also: gsp_igft, gsp_gwft, gsp_compute_fourier_basis
%

% Author : Nathanael Perraudin, David I Shuma
% Date: 

s = size(f);

if ~isnumeric(G);
    if ~gsp_check_fourier(G)
       error(['GSP_GFT: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    %U = G.U;
    f_hat=reshape(G.U'*reshape(f,G.N,[]),s);
else
    %U = G;
    f_hat=reshape(G'*reshape(f,G.N,[]),s);
end


