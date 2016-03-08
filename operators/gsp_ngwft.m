function [ C ] = gsp_ngwft(G,f,g, param)
%GSP_NGWFT Normalized graph windowed Fourier transform
%   Usage:  G = gsp_ngwft(G,f,g, param);
%           G = gsp_ngwft(G,f,g);
%
%   Input parameters:
%         G     : Graph
%         f     : Graph signal
%         g     : window
%         param : Structure of optional parameter
%   Output parameters:
%         C     : Coefficient
%
%   This function compute the normalized graph windowed Fourier transform
%   of a signal f with the window g. The function returns a matrix of
%   size N^2*N.  
%
%   param a Matlab structure containing the following fields:
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%      param.lowmemory : use less memory. By default, it is 1.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_ngwft.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.1
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
% Testing: test_ngwft

% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lowmemory'), param.lowmemory=1 ; end

if ~isfield(G,'U')
   error(['GSP_NGWFT: You need first to compute the Fourier basis\n',...
       'You can do it with the function gsp_compute_fourier_basis']);
end

if ~param.lowmemory
    % Compute the Frame into a big matrix
    Frame=gsp_ngwft_frame_matrix(G,g,param);

    C=Frame'*f;

    C=reshape(C,G.N,G.N);
else
    % Compute the translate of g
    ghat=G.U'*g;
    Ftrans=sqrt(G.N)*G.U*(repmat(ghat,1,G.N).*G.U');
    C=zeros(G.N);
    for ii=1:G.N
         atoms = (repmat(1./G.U(:,1),1,G.N) .* ...
             G.U.*repmat(Ftrans(:,ii),1,G.N))';
         % normalization
         atoms = atoms./ repmat(sqrt(sum(abs(atoms).^2,2)),1,G.N);
         C(:,ii)=atoms*f;
    end
    
end


