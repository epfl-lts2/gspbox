function [W,filtertype] = gsp_jtv_design_dgw(G,K,psi_graph,psi_time)
%GSP_JTV_DESIGN_DGW Design a dynamic graph wavelet with separable kernels
%   Usage: [W,filtertype] = gsp_jtv_design_dgw(G,K);
%          [W,filtertype] = gsp_jtv_design_dgw(G,K,psi,phi);
%
%   Input parameters:
%       G           : Time-Vertex Graph structure
%       K           : Time-Vertex Kernels
%       psi_graph   : Graph Kernel (default (@(x)1)
%       psi_time    : Time Kernel (default (@(t)1)
%   Output parameters:
%       W           : Dynamic Graph Wavelet
%       filtertype  : Filter domain ts
%
%   Design a dynamic graph wavelet with separable kernels
%   kron(psi_graph,transpose(psi_time)) and time-vertex kernel K(L,t)
%
%   Additional parameters
%   ---------------------
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1)
%
%
%   Url: https://epfl-lts2.github.io/gspbox-html/doc/filters/gsp_jtv_design_dgw.html

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
% Date : July 2016

if nargin<4
    psi_graph = {@(x)1};
end
if nargin<3
    psi_time = {@(t)1};
end

if ~iscell(K)
    K={K};
end

if ~iscell(psi_graph)
    psi_graph={psi_graph};
end

if ~iscell(psi_time)
    psi_time={psi_time};
end

Nk=numel(K);
Nt=numel(psi_time);
Ng=numel(psi_graph);

W = cell(Nt,Ng,Ng);

if (Nk*Ng*Nt)>1
    for nk=1:Nk
        for ng=1:Ng
            for nt=1:Nt
                W{nt,ng,nk}=gsp_jtv_design_dgw(G,K{nk},psi_graph{ng},psi_time{nt});
            end
        end
    end
    W=vec(W);
    return
end


W = @(x,t) K{1}(x,t).*psi_graph{1}(x).*psi_time{1}(t);

filtertype = 'ts';
