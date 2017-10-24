function coords = gsp_lle(G, dim, param)
%GSP_LLE Local Linear Embedding
%   Usage: coords = gsp_lle(G, dim);
%          coords = gsp_lle(G, dim, param);
%
%   Input parameters
%         G          : Graph
%         dim        : Dimensionality of the embedding
%         param      : Structure of optional parameters
%
%   Output parameters
%         coords     : Coordinates of the embedding
%
%   This function uses the weight matrix of a graph G, in order to compute
%   a *dim* -dimensional embedding (output coordinates). The algorithm used
%   is Locally Linear Embedding (LLE). Warning, this function might not
%   work if the graph is not connected.
%
%   *param* is a structure with optional parameters:
%
%   * *param.tol*    : Tolerance for the spectral gap (default 1e-6).
%   * *param.kernel* : The kernel used to create the graph weight matrix:
%     + 'exp'        : exponential kernel ($e^(frac{-d_{ij}}{sigma^2})$)
%     + '1/x'        : inverse of x kernel ($frac{1}{sigma+d_{ij}}$)
%     + '1/x^2'      : inverse of x^2 kernel ($frac{1}{(sigma+d_{ij})^2}$)
%     + 'resistance' : Resistance distance.
%   * *param.k*      : Max number of nearest neighbors. If not defined, the
%
%   number of neighbors varies from node to node since the algorithm
%   considers all the columns $j$ of the weight matrix that have non
%   zero values on the line $i$ as the neighbors of $i$.
%
%   References: saul2000introduction
%
%   See also: gsp_update_coordinates gsp_isomap gsp_laplacian_eigenmaps
%
%   Demo: gsp_demo_graph_embedding

% Authors : Dion O. E. Tzamarias
% Date    : 20/11/2015

%% TODO Fix bug param.k
%%

if nargin<3
    param = struct;
end


if ~isfield(param,'tol'), param.tol = 1e-6; end
if ~isfield(param,'kernel'), param.kernel = 'exp'; end


% Handling default parameters

W = zeros(G.N,G.N);

% if ~isfield(param,'k')
    [D] = gsp_weight2distance(G);
    [rows, cols ] = find(D);
    idx = cell(cols(end),1);
    for jj=1:cols(end)
        idx{jj} = rows(cols==jj).';
    end
% else
%     [D,idx] = gsp_weight2distance(G, param.kernel, param.k);
% end
D = D.^2; % D contains the squared distances
for ii=1:length(D)
    
    %   C = 0.5*(repmat(sum(D(idx{ii},idx{ii}),2),1,length(idx{ii}))/length(idx{ii})...
    %       + repmat(sum(D(idx{ii},idx{ii}),2).',length(idx{ii}),1)/length(idx{ii}) ...
    %       - D(idx{ii},idx{ii}) - repmat(sum(sum(D(idx{ii},idx{ii}),1)),...
    %       length(idx{ii}),length(idx{ii}))/length(idx{ii}).^2);
    C = 0.5*(bsxfun(@plus,sum(D(idx{ii},idx{ii}),2),sum(D(...
        idx{ii},idx{ii}),1))/length(idx{ii})-bsxfun(@plus,...
        D(idx{ii},idx{ii}),sum(sum(D(idx{ii},idx{ii}),1))...
        /length(idx{ii}).^2));
    
    if abs(det(C))< param.tol  %   C is singular
        C = C + (param.tol*trace(C)) / length(idx{ii})...
            * eye(size(C));
    end
%     if abs(det(C))< param.tol*1e-3
%         warning('C is singular')
%     end
    
    W(ii,idx{ii}) = C\ones(length(idx{ii}),1);
    W(ii,:) = W(ii,:)./sum(W(ii,:));
end

M = (eye(G.N,G.N) - W)' * (eye(G.N,G.N) - W);
M = sparse(M);

% Compute first dim eigenvectors of M
options.disp = 0;
options.isreal = 1;
options.issym = 1;
% only need bottom (no_dims + 1) eigenvectors
tol = param.tol;
[coords, e] = eigs(M, dim + 1, tol, options);
e = diag(e);

if nnz(e<param.tol) > 1
    disp('Multiple zero eigenvalues')
end

[~ , ind] = sort(e, 'ascend');
coords = coords(:,ind);
coords = coords(:,2:end);
end


%   The input parameter type refers to the type of the problem:
%   * 'weight_matrix' : Solve problem using the weight matrix of the graph G
%     Use the pairwise distance equation in order to reconstruct data from
%     neighbours.
%   * 'coords'       : Slove the problem using the coordinates. This is the
%     typical LLE aproach.



% switch lower(type)
%     case 'coords'
%         for ii=1:G.N
%
%             k_near_indx = knnsearch(G.coords,G.coords(ii,:),'k',n_neighbor+1);
%             k_nearest = G.coords(k_near_indx(2:end),:);
%
%             C = (repmat(G.coords(ii,:),n_neighbor,1) - k_nearest)* ...
%                 (repmat(G.coords(ii,:),n_neighbor,1) - k_nearest)' ;
%
%             if abs(det(C))<tol  % C is not singular
%                 C = C + ((tol*trace(C))^2) / n_neighbor * eye(size(C));
%             end
%             if abs(det(C))<tol
%                 warning('C is not singular')
%             end
%
%             W(ii,k_near_indx(2:end)) = C\ones(n_neighbor,1);
%             W(ii,:) = W(ii,:)./sum(W(ii,:));
%         end
%
%     case 'weight_matrix'
