function [Gs]=gsp_graph_multiresolution_old(G,num_levels,varargin)
%GSP_GRAPH_MULTIRESOLUTION  Compute a multiresolution of graphs
%   Usage:  [Gs]=gsp_graph_multiresolution(G,num_levels);
%           [Gs]=gsp_graph_multiresolution(G,num_levels,param);
%
%   Input parameters:
%         G                          : Graph structure.
%         num_levels                 : Number of times to downsample and coarsen the graph.
%   Output parameters:
%         Gs                         : Cell array, with each element containing a graph structure represent a reduced graph.
%   Additional parameters:
%         param.downsampling_method  : The graph downsampling method (default='largest_eigenvector')
%         param.reduction_method     : The graph reduction method (default='Kron')
%         param.sparsify             : To perform a spectral sparsification step immediately after the graph reduction (default=1)
%         param.sparsify_epsilon     : Parameter epsilon used in the spectral sparsification (default=min(10/sqrt(G.N),.3))
%         param.compute_full_eigen   : To also compute the graph Laplacian eigenvalues and eigenvectors for every graph in the multiresolution sequence (default=0)
%         
%   'gsp_graph_multiresolution(G,num_levels)' computes a multiresolution of 
%   graph by repeatedly downsampling and performing graph reduction. The
%   default downsampling method is the largest eigenvector method based on 
%   the polarity of the components of the eigenvector associated with the 
%   largest graph Laplacian eigenvalue. The default graph reduction method
%   is Kron reduction followed by a graph sparsification step.
%
%   See also:  
%
%   Demos:  
% 
%   References: 
%
%   AUTHORS : David I Shuman, Elle Weeks, Andre Archer, Stefan Faridani, Yan Jin.
%   TESTING: 
%   REFERENCE:

if nargin>2
    param=varargin{1};
else
    param=0;
end

if ~isfield(param,'reduction_method')
    reduction_method='kron'; 
else
    reduction_method=param.reduction_method;
end

if ~isfield(param,'downsampling_method')
    downsampling_method='largest_eigenvector';
else
    downsampling_method=param.downsampling_method;
end

if ~isfield(param,'sparsify')
    sparsify=1; 
else
    sparsify=param.sparsify;
end


if ~isfield(param,'sparsify_epsilon')
    sparsify_epsilon=min(10/sqrt(G.N),.3);
else
    sparsify_epsilon=param.sparsify_epsilon;
end

if ~isfield(param,'compute_full_eigen')
    compute_full_eigen=0; 
else
    compute_full_eigen=param.compute_full_eigen;
end


%set up cell for multiresolutions of graphs
Gs=cell(num_levels+1,1);
Gs{1}=G;
if compute_full_eigen
    if (~isfield(Gs{1},'U') || ~isfield(Gs{1},'e') )
        Gs{1}=gsp_compute_fourier_basis(Gs{1});
    end
else
    if ~isfield(Gs{1},'lmax')
        Gs{1}=gsp_estimate_lmax(Gs{1});
    end
end
Gs{1}.idx=(1:Gs{1}.N)';
Gs{1}.orig_idx=Gs{1}.idx;


for lev=1:num_levels
    
    % Graph downsamping: get indices to keep for the new lower resolution graph
    switch downsampling_method
        case 'largest_eigenvector'
            if compute_full_eigen
                largest_eigenvector = Gs{lev}.U(:,Gs{lev}.N);
            else
                [largest_eigenvector,~]=eigs(Gs{lev}.L,1); 
            end
             largest_eigenvector = largest_eigenvector * sign(largest_eigenvector(1));
            
            nonnegative_logicals=(largest_eigenvector >= 0);
            if sum(nonnegative_logicals) == 0
                error('Too many pyramid levels. Try fewer.');
            end
            keep_inds=find(nonnegative_logicals==1);
            
        % we can add other downsampling methods here
        
        otherwise
            error('Unknown graph downsampling method');
    end
    
   
    % Graph reduction: rewire the new lower resolution graph to form weighted adjacency and Laplacian matrices
    switch reduction_method
        case 'kron'
            % Kron reduction
            Gs{lev+1}.L=gsp_kron_reduce_old(Gs{lev}.L,keep_inds);
          
        % we can add other graph reduction methods here
        
        otherwise
            error('Unknown graph reduction method');
    end
    
    % Spectral sparsification
    if sparsify
        N=size(Gs{lev+1}.L,1);
        sparsify_epsilon=max(sparsify_epsilon,2/sqrt(N));
    %    gsp_reset_seed();
        [Gs{lev+1}.L,~,~] = gsp_graph_sparsify_old(Gs{lev+1}.L,sparsify_epsilon);
    end
    
    % Create the new graph from the reduced weighted adjacency matrix 
    new_W=diag(diag(Gs{lev+1}.L))-Gs{lev+1}.L;
    Gs{lev+1}=gsp_graph(new_W);
    Gs{lev+1}=gsp_copy_graph_attributes(Gs{lev},'unknown',Gs{lev+1});

    % Copy the coordinates of the subsampled vertices
    Gs{lev+1}.coords = Gs{lev}.coords(keep_inds,:);
    
    % Update indexing
    Gs{lev+1}.idx=keep_inds;
    Gs{lev+1}.orig_idx=Gs{lev}.orig_idx(keep_inds);
    
    % Compute full eigendecomposition of new graph, if desired
    if compute_full_eigen
        Gs{lev+1}=gsp_compute_fourier_basis(Gs{lev+1});
    else
        Gs{lev+1}=gsp_estimate_lmax(Gs{lev+1});
    end
     
end


end



  

