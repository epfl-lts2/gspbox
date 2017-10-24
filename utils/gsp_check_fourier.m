function bool = gsp_check_fourier(G)
%GSP_CHECK_FOURIER Check is the Fourier basis is computed
%   Usage:  bool = gsp_check_fourier(G):
%
%   Input parameters:
%       G           : Graph structure
%   Output parameters:
%       bool        : boolean
%
%   This function check is the Laplacian eigenvalues and eigenfunctions are
%   computed.
%

bool = ( isfield(G,'e') && isfield(G,'U') );

end