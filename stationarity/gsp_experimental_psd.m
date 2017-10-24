function gf = gsp_experimental_psd(G,C)
%GSP_EXPERIMENTAL_PSD Experimental power density function
%   Usage:  f_hat=gsp_gft(G,f);
%
%   Input parameters:
%         G          : Graph or Fourier basis
%         C          : Covariance matrix
%   Output parameters:
%         gf         : PSD filter
%
%   This function estimate the PSD from the covariance matrix with
%
%   .. T = U' * C * U 
%
%   .. math:: T = U^{*} C U
%
%   where $U$ is the graph Fourier basis. The function then interpolates
%   the diagonal of $T$ with splines.
%
%   To compute the Fourier basis of a graph G, you can use the function::
%
%           G = gsp_compute_fourier_basis(G);
%
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016


if ~isnumeric(G);
    if ~isfield(G,'U')
       error(['GSP_EXPERIMENTAL_PSD: You need first to compute the Fourier basis\n',...
           'You can do it with the function gsp_compute_fourier_basis']);
    end
    U = G.U;
else
    U = G;
end


    Cf = U'*C*U;
    
    cf = real(diag(Cf));
    [ e,ind ]= unique(G.e);
    gf = @(s) max(spline(e,cf(ind),s),0);
end

% % Old code with polynomial fitting
%     X = zeros(N,Order+1);
%     X(1,:) = 1;
%     for ii = 1:Order
%         X(ii+1) = X(ii+1) .* G.e;
%     end
%     
%     alpha = Pinv(X)*cf;
%     [p,S,mu] = polyfit(G.e,cf,Order);
%     gf = @(x) abs(polyval(p,x,S,mu))+ eps;

