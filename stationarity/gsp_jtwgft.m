function jtwgft = gsp_jtwgft(G, w, s, param)
%GSP_JTWGFT Compute the joint vertex-time windowed graph fourier transform
%   Usage:  jtwgft = gsp_jtwgft(G, w, s)
%           jtwgft = gsp_jtwgft(G, w, s, param)
%
%   Input parameters:
%         G          : Graph structure
%         w          : windows
%         s          : vertex-time signal
%         param      : structure of optional parameters
%   Output parameters:
%         jtwgft     : Joint vertex-time windowed graph Fourier transform
%
%   Additional parameters
%   ---------------------
%   * *param.boundary*    : 'periodic' or 'reflecting' (default 'periodic')
%   * *param.a*           : hop size in time (default 1)
%   * *param.M*           : number of frequency (default size(s,2))
%

% Author : Nathanael Perraudin
% Date   : 30 April 2016

if nargin<4
    param=struct;
end
if ~isfield(param,'boundary'), param.boundary = 'periodic'; end
if ~isfield(param,'a'), param.a = 1; end
if ~isfield(param,'M'), param.M = size(s,2); end


switch param.boundary
    case 'periodic'
        % 1) perform DGT
        S = dgt(transpose(s),w,param.a,param.M);
        % 2) Reshape
        S = reshape(S,[],G.N);
        % 3) GFT
        Shat = transpose(gsp_gft(G,transpose(S)));
        % 4) Reshape
        jtwgft = reshape(Shat,param.M,[],G.N);
    case 'reflecting'
        % 1) perform DGT
        S = wmdct(transpose(s),w,param.M,size(s,2));
        % 2) Reshape
        S = reshape(S,[],G.N);
        % 3) GFT
        Shat = transpose(gsp_gft(G,transpose(S)));
        % 4) Reshape
        jtwgft = reshape(Shat,param.M,[],G.N);
    otherwise
        error('Unknown boundary condition');
end

end