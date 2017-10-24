function c = gsp_jtv_frame_analysis(F,s)
%GSP_JTV_FRAME_ANALYSIS Analysis operator in time-vertex domain
%   Usage:  c = gsp_jtv_frame_analysis(F,s,param);
%
%   Input parameters:
%         F          : Frame matrix (vertex_loc x time_loc x vertex x time)
%         s          : Time-Vertex signal
%   Output parameters:
%         c          : Coefficients matrix
%   ---------------------
%   This function compute the analysis operator using frame matrix

% Author :  Francesco Grassi
% Date   : July 2016

[N,lag,~,T]=size(F);
Ts = size(s,2);

if Ts<(T+1)/2
    s = [ zeros(N,T-Ts-1) s ];
elseif T<Ts
    error('Atoms time length should be bigger or equal to the signal time length')
end


c = reshape(reshape(F,N*lag,[])*s(:),N,[]);