function [ ft ] = gsp_translate_old(G, f, i)
%GSP_TRANSLATE Generalized translation of the signal f to the node i
%   Usage: ft = gsp_translate(G, f, i);
%
%   Input parameters
%       G   : Graph
%       f   : Signal (column)
%       i   : Indices of vertex (int)
%   Output parameters
%       ft  : translate signal
%
%   This function translate the column vector f onto the node i. If f*
%   is a matrix, the translation will be done to each column.
%
%

fhat=gsp_gft(G,f);
nt = size(f,2);

ft = sqrt(G.N)*gsp_igft(G,fhat .* ...
    repmat(transpose(G.U(i,:)),1,nt));


end