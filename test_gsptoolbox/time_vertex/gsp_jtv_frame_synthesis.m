function s = gsp_jtv_frame_synthesis(F,c)
%GSP_JTV_FRAME_SYNTHESIS Synthesis operator in time-vertex domain
%   Usage:  s = gsp_jtv_frame_synthesis(F,c,param);
%
%   Input parameters:
%         F          : Frame matrix (vertex_loc x time_loc x vertex x time)
%         c          : Coefficients matrix
%   Output parameters:
%         s          : Time-Vertex signal
%         
%   This function compute the synthesis operator using frame matrix

% Author :  Francesco Grassi
% Date   : July 2016


[N,lag,~,T]=size(F);
s = reshape(reshape(F,N*lag,[])'*c(:),N,[]);