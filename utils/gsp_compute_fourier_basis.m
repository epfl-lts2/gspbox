function [G] = gsp_compute_fourier_basis(G,param)
%GSP_COMPUTE_FOURIER_BASIS Compute the fourier basis of the graph G
%   Usage:  G = gsp_compute_fourier_basis(G);
%           G = gsp_compute_fourier_basis(G,param);
%
%   Input parameters:
%         G          : Graph structure (or cell array of graph structure) 
%         param      : structure of optional parameters
%   Output parameters:
%         G          : Graph structure (or cell array of graph structure)
%
%   'gsp_compute_fourier_basis(G)' computes a full eigendecomposition of the graph
%   Laplacian G.L:
%
%      L = U Lambda U* 
%
%   where Lambda is a diagonal matrix of the Laplacian eigenvalues. 
%   G.e is a column vector of length G.N containing the Laplacian
%   eigenvalues. The function will store the basis U, the eigenvalues
%   e, the maximum eigenvalue lmax and G.mu the coherence of the
%   Fourier basis into the structure G.
% 
%   Example:
%
%       N = 50;
%       G = gsp_sensor(N);
%       G = gsp_compute_fourier_basis(G);
%       gsp_plot_signal(G,G.U(:,2));
% 
%   References:
%     F. R. K. Chung. Spectral Graph Theory. Vol. 92 of the CBMS Regional
%     Conference Series in Mathematics, American Mathematical Society, 1997.
%     
%     
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_compute_fourier_basis.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

% Author : David I Shuman, Nathanael Perraudin
% Testing: test_operators

if nargin < 2
    param = struct;
end


if numel(G)>1
    Ng = numel(G);
    for ii = 1:Ng
       G{ii} = gsp_compute_fourier_basis(G{ii}, param);
    end     
    return;
end

if ~isfield(param,'verbose'), param.verbose = 1; end




if gsp_check_fourier(G)
    if param.verbose
        warning(['Laplacian eigenvalues or eigenvectors ',...
            'are already associated with this graph']);
    end
end

if G.N > 15000
    if param.verbose
        error('Too big matrix to perform full eigenvalue decomposition.'); 
    end
end

if G.N > 3000
    if param.verbose
        warning(['Performing full eigendecomposition ',...
            'of a large matrix may take some time...']); 
    end
end
    
if isfield(G,'type') &&  strcmp(G.type,'ring')==1 % && mod(G.N,2)==0 
    U = dftmtx(G.N)/sqrt(G.N);
    E = (2-2*cos(2*pi*(0:G.N-1)'/G.N));
    inds = gsp_classic2graph_eig_order( G.N );
%     [G.E, inds]=sort(E,'ascend');
    G.e = E(inds);
    if strcmp(G.lap_type,'normalized')
        G.e = G.e/2;
    end
    
    G.U = U(:,inds);
else
    if ~isfield(G,'L')
        error('Graph Laplacian is not provided.');
    end
    [G.U, G.e] = gsp_full_eigen(G.L);
end

G.lmax=max(G.e);

if isfield(G,'Gm')
    G = gsp_compute_oose_fourier_basis(G);
end



G.mu = max(abs(G.U(:)));

end


function [U,E] = gsp_full_eigen(L)
%GSP_FULL_EIGEN Compute and order the eigen decomposition of L

    % Compute and all eigenvalues and eigenvectors 
%     try
%         [eigenvectors,eigenvalues]=eig(full(L+L')/2);
%     catch
        [eigenvectors,eigenvalues,~]=svd(full(L+L')/2);
%     end
    
    % Sort eigenvectors and eigenvalues
    [E,inds] = sort(diag(eigenvalues),'ascend');
    eigenvectors=eigenvectors(:,inds);
    
    % Set first component of each eigenvector to be nonnegative
    signs=sign(eigenvectors(1,:));
    signs(signs==0)=1;
    U = eigenvectors*diag(signs);
end

function D = dftmtx(n)

n = signal.internal.sigcasttofloat(n,'double','dftmtx','N',...
  'allownumeric');

D = fft(eye(n));

end

