function r = gsp_stationarity_ratio(G, C, param)
%GSP_STATIONARITY_RATIO This is a number from [0 to 1] depicting how close to the PCA basis the graph basis is
%
%   The method examines the percentage of the data variance that is not in
%   the diagonal of the covariance of the GFT of X. The index can be used to
%   describe how well a graph fits a given data matrix X, or distribution. 
%   An index of 0 means that the data are graph stationary on G.
%
%   Usage:  r = gsp_stationarity_ratio(G, C)
%
%   Input parameters:
%         G          : Graph
%         C          : Covariance matrix
%   Output parameters:
%         r          : Ratio
%
%   Optional parameters: 
%   params.verbose   :0 = nothing, 1 = plot the covariance of GFT(x) (default 0)
% 
%   This function compute the ratio of energy contained into the diagonal
%   of the Fourier covariance matrix:
%
%   .. T = U' * C * U 
%
%   .. math:: T = U^{*} C U
%
%   References: perraudin2016stationary


% Authors : Nathanael Perraudin, Andreas Loukas
% Date    : 6 January 2016

if nargin < 3, param = struct(); end

if ~isfield(param,'verbose'), param.verbose = 0; end

if not(isfield(G, 'U'))
    G = gsp_compute_fourier_basis(G);
end


CF = G.U' * C * G.U;

r = gsp_diagonal_ratio(CF);


if param.verbose > 0
    figure; imagesc(abs(CF)); 
    title('G.U^*  * (X * X^*) * G.U')
end

end


function r = gsp_diagonal_ratio(M)

if size(M,1) == 1 || size(M,2) ==1
    error('This function acts on martrices.')
end

r = norm(diag(M))/norm(M,'fro');

end