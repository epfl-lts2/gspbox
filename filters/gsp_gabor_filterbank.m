function [ g ] = gsp_gabor_filterbank( G,k )
%GSP_GABOR_FILTERBANK returns a Gabor filterbank
%   Usage: g = gsp_gabor_filterbank( G,k );
%
%   Input parameters:
%       G       : Graph
%       k       : kernel
%
%   Output parameters:
%       g       : filterbank
%   
%   This function creates a filterbank with the kernel *k*. Every filter is
%   centered in a different frequency
%

% Author: Nathanael Perraudin
% Date  : 13 June 2014

error('To rename and changed!')
if ~isfield(G,'e')
   error(['GSP_GABOR_FILTERBANK: You need first to compute the Fourier basis\n',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

Nf = length(G.e);

g = cell(Nf,1);

for ii=1:Nf
    g{ii} = @(x) k(x-G.e(ii));
end


end

