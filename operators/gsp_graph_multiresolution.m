function [Gs]=gsp_graph_multiresolution(G,num_levels,param)
%GSP_GRAPH_MULTIRESOLUTION  Compute a multiresolution of graphs
%   Usage:  [Gs]=gsp_graph_multiresolution(G,num_levels);
%           [Gs]=gsp_graph_multiresolution(G,num_levels,param);
%
%   Input parameters:
%         G           : Graph structure
%         num_levels  : Number of times to downsample and coarsen the graph
%         param       : Optional structure of parameters
%   Output parameters:
%         Gs          : Cell array of graphs
%         
%   'gsp_graph_multiresolution(G,num_levels)' computes a multiresolution of 
%   graph by repeatedly downsampling and performing graph reduction. The
%   default downsampling method is the largest eigenvector method based on 
%   the polarity of the components of the eigenvector associated with the 
%   largest graph Laplacian eigenvalue. The default graph reduction method
%   is Kron reduction followed by a graph sparsification step.
%   *param* is a structure of optional parameters containing the following
%   fields:
%
%   * *sparsify*: To perform a spectral sparsification step immediately
%     after the graph reduction (default=1) 
%   * *sparsify_epsilon*: Parameter epsilon used in the spectral
%     sparsification (default=min(10/sqrt(G.N),.3))   
%   * *downsampling_method*: The graph downsampling method
%     (default='largest_eigenvector') 
%   * *reduction_method*: The graph reduction method (default='Kron')
%   * *compute_full_eigen*: To also compute the graph Laplacian eigenvalues
%     and eigenvectors for every graph in the multiresolution sequence
%     (default=0)  
%
%   Example:::
%
%             N = 500;
%             G = gsp_sensor(N);
%             Nlevel = 5;
% 
%             Gs = gsp_graph_multiresolution(G, Nlevel);
% 
%             figure;
%             for ii = 1:numel(Gs)
%                 subplot(2,3,ii)
%                 gsp_plot_graph(Gs{ii})
%                 title(['Reduction level: ', num2str(ii-1)]);
%             end
%
%   See also: gsp_pyramid_analysis gsp_pyramid_synthesis gsp_pyramid_cell2coeff
%
%   Demo: gsp_demo_pyramid
%
%   References: shuman2013framework
%
%   Authors : David I Shuman, Elle Weeks, Andre Archer, Stefan Faridani, Yan Jin, Nathanael Perraudin.
%   Date: 26 November 2015
%   Testing: test_operators 


if nargin < 3
    param = struct;
end

if ~isfield(param,'sparsify'), param.sparsify = 1; end;
if ~isfield(param,'compute_full_eigen'), param.compute_full_eigen = 0; end;
if ~isfield(param,'sparsify_epsilon'), param.sparsify_epsilon = min(10/sqrt(G.N),.3); end;

if ~isfield(param,'reduction_method')
    param.reduction_method='kron'; 
end

if ~isfield(param,'downsampling_method')
    param.downsampling_method='largest_eigenvector';
end


if param.compute_full_eigen
    if (~isfield(G,'U') || ~isfield(G,'e') )
        G=gsp_compute_fourier_basis(G);
    end
else
    if ~isfield(G,'lmax')
        G=gsp_estimate_lmax(G);
    end
end

%set up cell for multiresolutions of graphs
Gs=cell(num_levels+1,1);
Gs{1}=G;
Gs{1}.mr.idx=(1:Gs{1}.N)';
Gs{1}.mr.orig_idx=Gs{1}.mr.idx;

if param.compute_full_eigen
    if (~isfield(Gs{1},'U') || ~isfield(Gs{1},'e') )
        Gs{1}=gsp_compute_fourier_basis(Gs{1});
    end
else
    if ~isfield(Gs{1},'lmax')
        Gs{1}=gsp_estimate_lmax(Gs{1});
    end
end

for lev=1:num_levels
    Gs{lev+1}.directed=0;
    
    % Graph downsampling: get indices to keep for the new lower resolution graph
    switch param.downsampling_method
        case 'largest_eigenvector'
            if isfield(Gs{lev},'U')
                largest_eigenvector = Gs{lev}.U(:,Gs{lev}.N);
            else
                [largest_eigenvector,~]=eigs(Gs{lev}.L,1); 
            end
            largest_eigenvector=largest_eigenvector*sign(largest_eigenvector(1));
            nonnegative_logicals=(largest_eigenvector >= 0);
            if sum(nonnegative_logicals) == 0
                error('Too many pyramid levels. Try fewer.');
            end
            keep_inds=find(nonnegative_logicals);
            
        % we can add other downsampling methods here
        
        otherwise
            error('Unknown graph downsampling method');
    end
    
   
    % Graph reduction: rewire the new lower resolution graph to form weighted adjacency and Laplacian matrices
    switch param.reduction_method
        case 'kron'
            % Kron reduction
            Gs{lev+1}.L=gsp_kron_reduce(Gs{lev}.L,keep_inds);
          
        % we can add other graph reduction methods here
        
        otherwise
            error('Unknown graph reduction method');
    end
 
    % Create the new graph from the reduced weighted adjacency matrix 
    Gs{lev+1}.N=size(Gs{lev+1}.L,1);
    Gs{lev+1}.W=diag(diag(Gs{lev+1}.L))-Gs{lev+1}.L;
    Gs{lev+1}.A=sign(Gs{lev+1}.W);
    Gs{lev+1}=gsp_copy_graph_attributes(Gs{lev},1,Gs{lev+1});
    
    % Spectral sparsification
    if param.sparsify
        if Gs{lev+1}.N>2
            sparsify_epsilon=max(param.sparsify_epsilon,2/sqrt(Gs{lev+1}.N));
            Gs{lev+1} = gsp_graph_sparsify(Gs{lev+1},sparsify_epsilon);
        end
    end
    
    % Copy the coordinates of the subsampled vertices
    Gs{lev+1}.coords = Gs{lev}.coords(keep_inds,:);
    Gs{lev+1}.type='from multiresolution';
        
    % Update indexing
    Gs{lev+1}.mr.idx=keep_inds;
    Gs{lev+1}.mr.orig_idx=Gs{lev}.mr.orig_idx(keep_inds);
    
    % Compute full eigendecomposition of new graph, if desired
    if param.compute_full_eigen
        Gs{lev+1}=gsp_compute_fourier_basis(Gs{lev+1});
    else
        Gs{lev+1}=gsp_estimate_lmax(Gs{lev+1});
    end
     
end


end



  

