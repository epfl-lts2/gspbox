function [f, label] = gsp_jtv_fa(G,shift)
%GSP_JTV_FA Frequency axis for the joint time vertex framework
%   Usage:  f = gsp_jtv_fa(G);
%
%   Input parameters:
%       G       : Time-vertex graph structure
%       shift   : Boolean value: 1 to apply fftshift on the frequency vector, 0 otherwise
%   Ouput parameters:
%       f       : Frequency axis (row vector)
%

% Author: Nathanael Perraudin, Francesco Grassi
% Date  : September 2016

if nargin<2
    shift = 0;
end

if ~gsp_check_jtv(G)
    error('GSP_JTV_FA needs the time dimension. Use GSP_JTV_GRAPH');
end

if isempty(G.jtv.NFFT)
    NFFT = G.jtv.T;
else
    NFFT = G.jtv.NFFT;
end

if shift
    f = fftshift(gsp_cfa( NFFT,G.jtv.fs )');
else
    f = gsp_cfa( NFFT,G.jtv.fs )';
end




label = '\omega';

end