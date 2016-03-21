function [ G ] = gsp_create_laplacian( G,type )
%GSP_CREATE_LAPLACIAN create the graph laplacian of the graph G
%   Usage: G = gsp_create_laplacian( G,type );
%          G = gsp_create_laplacian( G );
%
%   Input parameters:
%       G   : Graph structure (or cell array of graph structure) 
%       type: Type of laplacian (string)
%   Output parameters:
%       G   : Graph structure (or cell array of graph structure) 
%
%   This function create the graph laplacian of the graph G and store it
%   into G.
%
%   The variable type contains the different laplacian type. For
%   undirected graph, the following type are availlable:
%
%    combinatorial*: Non normalized laplacian. This is the default.
%
%          L =  D  - W 
%
%   And for directed graph, the following types are availlable.
%
%    combinatorial : Non normalized laplacian. This is the default
%
%          L =  1/2 [ D^+ + D^- - W - W^*]
%
%    chung*: Normalized laplacian with the Perron eigenvector
%
%        L_cn = I - 1/2 [Pi^0.5 P Pi^-0.5 + Pi^-0.5 P^T Pi^0.5 ]
%
%       
%   References:
%     F. Chung. Laplacians and the cheeger inequality for directed graphs.
%     Annals of Combinatorics, 9(1):1--19, 2005.
%     
%     
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_create_laplacian.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.5.2
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


% Author: Nathanael Perraudin
% Date  : 09.12.2013

if numel(G)>1
    Ng = numel(G);
    for ii = 1:Ng
        if nargin<2
            G{ii} = gsp_create_laplacian(G{ii});
        else
            G{ii} = gsp_create_laplacian(G{ii}, type);
        end
    end     
    return;
end


if nargin<2
    if ~isfield(G,'lap_type')
        type='combinatorial';
        G.lap_type = type;
    else
        type = G.lap_type;
    end
end

if isfield(G,'hypergraph') && G.hypergraph
    G.de = sum(G.W >0,1)';
    G.dv = sum(G.W.^2,2);
    
%     switch type
%         case 'normalized'        
            G.A = G.W*G.W' - diag(G.dv);
            G.L = eye(G.N) - ...
                diag(G.dv.^(-0.5)) * G.W * ...
                diag(G.de.^(-1)) * ...
                G.W' *diag(G.dv.^(-0.5));
            G.lap_type = 'normalized';
%         otherwise
%             error('Unknown laplacian type')
%     end
    return
    
end

if G.directed
    D1 = sum(G.W,2);
    D2 = sum(G.W,1);
    
    switch type
        case 'combinatorial'
            G.L=0.5 * (diag(D1) + diag(D2) - G.W - G.W');
        case 'chung'
            [phi,P] = compute_perron(G.W);
            Phiup=diag(sparse(phi.^(0.5)));            
            Phidw=diag(sparse(phi.^(-0.5)));       
            G.L = speye(G.N) - 0.5 * (Phiup * P  * Phidw + Phidw * P'  * Phiup );
            % Save the results in G
            G.P=P;
            G.phi=phi;
        case 'chung-non-normalized'
            [phi,P] = compute_perron(G.W);
            Phiup=diag(sparse(phi.^(0.5)));            
            Phidw=diag(sparse(phi.^(-0.5))); 
            G.L = Phiup*(sparse(eye(G.N)) - 0.5 * (Phiup * P  * Phidw + Phidw * P'  * Phiup ))*Phiup;
            % Save the results in G
            G.P=P;
            G.phi=phi;
        case 'normalized'
            error('Not implemented yet. Ask Nathanael')
        case 'none'
            G.L=sparse(0);
        otherwise
            error(' Unknown laplacian type')
    end    
    % To avoid numerical errors of symetry
    G.L = (G.L +G.L')/2;
else
    D = sum(G.W,2);
    switch type
        case 'combinatorial'
            G.L=diag(D)-G.W;
        case 'normalized'
            Dn = diag(D.^(-0.5));
            G.L=speye(G.N)-Dn*G.W*Dn;
        case 'none'
            G.L=sparse(0);
        otherwise
            error(' Unknown laplacian type')
    end
end


if isfield(G,'Gm')
    G.Gm = gsp_create_laplacian(G.Gm,type);
end

% Update problematic fields
if isfield(G,'U')
    G = gsp_compute_fourier_basis(G);
end

if isfield(G,'Diff')
    G = rmfield(G,'Diff');
    G = gsp_adj2vec(G);
end

G.lap_type = type;


end


function [phi,P] = compute_perron(A)

%     N = size(A,1);

    % Remove the diagonal
    A=A-diag(diag(A));

    % Compute the Probablility matrix
    %P=A./repmat(sum(A,2),1,N);
%     P = bsxfun(@rdivide,A,sum(A,2));
    P = diag(sparse(1./sum(A,2)))*A;

    % Compute the perron vector of P
    [phi,max_eig_P] = eigs(P',1);
    % test if max_eig_P==1
    if abs(max_eig_P-1)>10e3*eps;
        fprintf(['\n  ---  Warning! The maximum eigenvalue of the probability ' ...
            'matrix is not 1 but %f  --- \n'],max_eig_P);
    end
    % Test if the perron vector is positive
    if sum(phi)<0; 
        phi=-phi;
    end
    
    
%     phi(phi<0) = -phi(phi<0);
%     warning('Things to be done here!!')
    if sum(phi<=10e3*eps)
        fprintf(['\n  ---  Warning! The perron vector has negative or '...
          'null entrie(s).\n       Is the graph strongly connected?  ---\n']);

    end
    % Normalization of phi
    phi=phi/norm(phi,1);

end

