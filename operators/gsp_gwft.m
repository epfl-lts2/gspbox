function [ C ] = gsp_gwft(G,g,f, param )
%GSP_GWFT Graph windowed Fourier transform
%   Usage:  C = gsp_gwft(G,g,f, param );
%           C = gsp_gwft(G,g,f);
%
%   Input parameters:
%         G     : Graph
%         g     : Window (graph signal or kernel)
%         f     : Graph signal (column vector)
%         param : Structure of optional parameter
%   Output parameters:
%         C     : Coefficient.
%
%   This function compute the graph windowed Fourier transform of a signal
%   f with the window g. The function returns a matrix of size N^2*N. 
%
%   param a Matlab structure containing the following fields:
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%      param.lowmemory : use less memory. By default, it is 1.
%
%   Reference: shuman2013windowed
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_gwft.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.0
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

% Author : Nathanael Perraudin
% Testing: test_gwft


% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lowmemory'), param.lowmemory=1 ; end

if ~isfield(G,'U')
   error(['GSP_GWFT: You need first to compute the Fourier basis\n',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

Nf = size(f,2);

if iscell(g)
    g = gsp_igft(G,g{1}(G.e));
end


if isa(g, 'function_handle')
    g = gsp_igft(G,g(G.e));
end

if ~param.lowmemory
    % Compute the Frame into a big matrix
    Frame=gsp_gwft_frame_matrix(G,g,param);

    C=Frame'*f;

    C=reshape(C,G.N,G.N,Nf);
else
    % Compute the translate of g
    ghat=G.U'*g;
    Ftrans=sqrt(G.N)*G.U*(repmat(ghat,1,G.N).*G.U');
    C=zeros(G.N);
    for jj = 1:Nf
        for ii=1:G.N
             C(:,ii,jj)=(repmat(1./G.U(:,1),1,G.N) .* ...
                      G.U.*repmat(Ftrans(:,ii),1,G.N))'*f(:,jj);
        end
    end
    
end


