function G = gsp_create_laplacian(G, lap_type)
%GSP_CREATE_LAPLACIAN create the graph laplacian of the graph G
%   Usage: G = gsp_create_laplacian(G, type);
%          G = gsp_create_laplacian(G);
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
%   The variable *type* contains the different laplacian type. For
%   undirected graph, the following type are availlable:
%
%   * *combinatorial*: Non normalized laplacian. This is the default.
%
%     ..   L =  D  - W 
%
%     .. math:: L= D - W
%   * *normalized*: Normalized laplacian
%     .. L_n = I - D^-0.5 W D^-0.5 
%
%     .. math:: L_{n}=I - D^{-\frac{1}{2}}WD^{-\frac{1}{2}}
%   * *none*: No laplacian
%
%   And for directed graph, the following types are availlable.
%
%   * *combinatorial* : Non normalized laplacian. This is the default
%
%     ..   L =  1/2 [ D^+ + D^- - W - W^*]
%
%     .. math:: L=\frac{1}{2} ( D^+ + D^- - W - W^*)
%
%   * *chung*: Normalized laplacian with the Perron eigenvector
%
%     .. L_cn = I - 1/2 [Pi^0.5 P Pi^-0.5 + Pi^-0.5 P^T Pi^0.5 ]
%
%     .. math:: L_{cn}=I-\frac{\Pi^{\frac{1}{2}}P\Pi^{-\frac{1}{2}}+\Pi^{-\frac{1}{2}}P^{*}\Pi^{\frac{1}{2}}}{2}
%   * *chung-non-normalized*: Experimental version of non normalized Chung
%     laplacian
%     .. L_c = Pi^0.5 L_cn Pi^0.5 
%
%     .. math:: L_{c}= \Pi^{\frac{1}{2}} L_{cn} \Pi^{\frac{1}{2}}
%
%
%   see also: gsp_laplacian
%       
%   References: chung2005laplacians


% Author: Nathanael Perraudin
% Date  : 09.12.2013, revised jun 2016

if numel(G)>1
    Ng = numel(G);
    for ii = 1:Ng
        if nargin<2
            G{ii} = gsp_create_laplacian(G{ii});
        else
            G{ii} = gsp_create_laplacian(G{ii}, lap_type);
        end
    end     
    return;
end


if nargin<2
    if ~isfield(G,'lap_type')
        lap_type='combinatorial';
        G.lap_type = lap_type;
    else
        lap_type = G.lap_type;
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
    
    switch lap_type
        case 'combinatorial'
            G.L=(diag(D1) + diag(D2) - G.W - G.W');
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
    % To avoid numerical errors of symmetry
    G.L = (G.L +G.L')/2;
else
    G.L = gsp_laplacian(G.W, lap_type);
end


if isfield(G,'Gm')
    G.Gm = gsp_create_laplacian(G.Gm,lap_type);
end

% Update problematic fields
if gsp_check_fourier(G)
    disp('Recomputing the Fourier matrix')
    G = rmfield(G,'U');
    G = rmfield(G,'e');
    G = rmfield(G,'lmax');
    G = gsp_compute_fourier_basis(G);
elseif isfield(G,'lmax')
    disp('Recomputing the maximum eigenvalue')
    G = gsp_estimate_lmax(G);
end



if isfield(G,'Diff')
    G = rmfield(G,'Diff');
    G = gsp_adj2vec(G);
end

G.lap_type = lap_type;


end


function [phi,P] = compute_perron(A)

%     N = size(A,1);

    % Remove the diagonal
    A = zero_diag(A);

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
