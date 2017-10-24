function F = gsp_jtv_compute_frame(G,g,filtertype)
% GSP_JTV_COMPUTE_FRAME Compute the frame tensor of the kernel g on the time-vertex graph G
%   Usage:  F = gsp_jtv_compute_frame(G,g);
%
%   Input parameters:
%         G          : Time-Vertex Graph structure
%         g          : Cell array of time-vertex filters
%         filtertype : Filter domain (ts,js,ts-array,js-array)
%   Output parameters:
%         F          : Frame tensor (vertex_loc x time_loc x vertex x time)
%
%   Compute the frame tensor of the kernel g on the time-vertex graph G
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_compute_frame.html

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this toolbox please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
% http://arxiv.org/abs/1408.5781

% Author :  Francesco Grassi
% Date : September 2016


if or(~isfield(G.jtv,'T'),~isfield(G.jtv,'fs'))
    error('GSP_JTV_COMPUTE_FRAME need the time dimension. Use GSP_JTV_GRAPH.');
end

if ~iscell(g)
    g = {g};
end

Nf = numel(g);

T   = G.jtv.T;
lag = G.jtv.lag;
fs  = G.jtv.fs;
N   = G.N;


F = zeros(N,lag,N,T,Nf);
param.lag = G.jtv.extension;

for n = 1:Nf
    for ii=1:N
        for jj=1:lag
            F(ii,jj,:,:,n) = gsp_jtv_filter_synthesis(G,g{n},filtertype,gsp_jtv_delta(G,ii,jj-T*param.lag,param));
        end
    end
    
end


end




%OLD CODE
% t = 0:1/fs:(T-1)/fs;
% tau = -(T-1):1/fs:(T-1);
%
% for n = 1:Nf
%     n
%     for jj=1:2*T-1
%         for ii = 1:T
%
%             tt=t(ii)-tau(jj);
%
%             h = gsp_jtv_filterbank(g,tt);
%
%             psi = real(gsp_filter_analysis(G,h{n},eye(N)));
%
%             F(:,jj,:,ii,n) = psi;
%
%         end
%     end
% end


