function [ fm ] = gsp_modulate( G,f,k )
%GSP_MODULATE Generalized modulation of the signal f to the frequency k
%   Usage: fm = gsp_modulate( G,f,k );
%
%   Input parameters
%       G   : Graph
%       f   : Signal (column)
%       k   : Indices of frequencies (int)
%   Output parameters
%       fm  : Modulated signal
%
%   This function modulate the column vector *f* onto the node i. If f is a
%   matrix, the modulation will be applicated to each column.
%

% Author: Nathanael Perraudin
% Date  : 09.12.2013

nt = size(f,2);

fm = sqrt(G.N)*repmat(f,1,nt).*repmat(G.U(:,k+1),1,nt);


end
