function [Gs,subsampled_vertex_indices]=gsp_tree_multiresolution(G,Nlevel,param)
%GSP_TREE_MULTIRESOLUTION  Compute a multiresolution of trees
%   Usage:  [Gs,subsampled_vertex_indices]=gsp_tree_multiresolution(G,num_levels);
%           [Gs,subsampled_vertex_indices]=gsp_tree_multiresolution(G,num_levels,param);
%
%   Input parameters:
%         G                          : Graph structure of a tree.
%         Nlevel                     : Number of times to downsample and coarsen the tree.
%   Output parameters:
%         Gs                         : Cell array, with each element containing a graph structure represent a reduced tree.
%         subsampled_vertex_indices  : Indices of the vertices of the previous tree that are kept for the subsequent tree.
%   Additional parameters:
%         param.root                 : The index of the root of the tree (default=1)
%         param.reduction_method     : The graph reduction method (default='resistance_distance')
%         param.compute_full_eigen   : To also compute the graph Laplacian eigenvalues for every tree in the sequence
%         
%   'gsp_tree_multiresolution(G,num_levels)' computes a multiresolution of 
%   trees by repeatedly downsampling and performing a graph reduction. The
%   downsampling is performed by keeping all vertices at even depths of the
%   tree from the root vertex. Options for the graph reduction method
%   include: 'unweighted', 'sum' (add the weight connecting a child
%   node to its parent and the weight connecting the parent to the
%   grandparent and use that weight for the edge connecting the child to
%   the grandparent in the new graph), or 'resistance_distance', which
%   preserves the resistance distances by setting the new weights according
%   to:
% 
%      W_ik =              1                      
%                     -----------
%                      1       1
%                     ---  +  ---  
%                    W_ij     W_jk
%
%   where W_{i,j} is the weight connecting a child to its parent in the
%   original tree, and W_{j,k} is the weight connecting the parent to the
%   grandparent in the original tree.
%
%   See also: gsp_kron_pyramid  
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/operators/gsp_tree_multiresolution.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.1
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

%   AUTHOR : David I Shuman, Nathanael Perraudin
%   TESTING: 
%   REFERENCE:

% TODO: transform this in tree pyramid

if nargin<3
    param = struct;
end


if ~isfield(param,'root'), 
    if isfield(G,'root')
        param.root = G.root;
    end
    param.root=1; end
if ~isfield(param,'reduction_method'), param.reduction_method='resistance_distance'; end
if ~isfield(param,'compute_full_eigen'), param.compute_full_eigen=0; end


root=param.root;
G.root = root;


Gs = cell(Nlevel+1,1);
Gs{1} = G;
if param.compute_full_eigen
    Gs{1} = gsp_compute_fourier_basis(Gs{1});
end
subsampled_vertex_indices = cell(Nlevel,1);

[depths,parents] = gsp_tree_depths(G.A,root);
old_W = G.W;


for lev=1:Nlevel
    
    % Identify the vertices in the even depths of the current tree
    down_odd = mod(round(depths),2);
    down_even = ones(Gs{lev}.N,1) - down_odd;
    keep_inds = find(down_even==1);
    subsampled_vertex_indices{lev} = keep_inds;

    % There will be one undirected edge in the new graph connecting each 
    % non-root subsampled vertex to its new parent. Here, we find the new
    % indices of the new parents
    [non_root_keep_inds,new_non_root_inds]=setdiff(keep_inds,root); 
    old_parents_of_non_root_keep_inds=parents(non_root_keep_inds);
    old_grandparents_of_non_root_keep_inds=parents(old_parents_of_non_root_keep_inds);
    new_non_root_parents=dsearchn(keep_inds,old_grandparents_of_non_root_keep_inds);
    
    % Create new weighted adjacency matrix via graph reduction
    [old_W_i_inds,old_W_j_inds,old_W_weights]=find(old_W);
    i_inds=[new_non_root_inds;new_non_root_parents];
    j_inds=[new_non_root_parents;new_non_root_inds];
    new_N=sum(down_even);
    switch param.reduction_method
        case 'unweighted'
            new_weights=ones(length(i_inds),1);
        case 'sum'
            old_weights_to_parents_inds=dsearchn([old_W_i_inds,old_W_j_inds],[non_root_keep_inds,old_parents_of_non_root_keep_inds]);
            old_weights_to_parents=old_W_weights(old_weights_to_parents_inds); %old_W(non_root_keep_inds,old_parents_of_non_root_keep_inds);
            old_weights_parents_to_grandparents_inds=dsearchn([old_W_i_inds,old_W_j_inds],[old_parents_of_non_root_keep_inds,old_grandparents_of_non_root_keep_inds]);
            old_weights_parents_to_grandparents=old_W_weights(old_weights_parents_to_grandparents_inds); %old_W(old_parents_of_non_root_keep_inds,old_grandparents_of_non_root_keep_inds);
            new_weights=old_weights_to_parents + old_weights_parents_to_grandparents;
            new_weights=[new_weights;new_weights];
        case 'resistance_distance'
            old_weights_to_parents_inds=dsearchn([old_W_i_inds,old_W_j_inds],[non_root_keep_inds,old_parents_of_non_root_keep_inds]);
            old_weights_to_parents=old_W_weights(old_weights_to_parents_inds); %old_W(non_root_keep_inds,old_parents_of_non_root_keep_inds);
            old_weights_parents_to_grandparents_inds=dsearchn([old_W_i_inds,old_W_j_inds],[old_parents_of_non_root_keep_inds,old_grandparents_of_non_root_keep_inds]);
            old_weights_parents_to_grandparents=old_W_weights(old_weights_parents_to_grandparents_inds); %old_W(old_parents_of_non_root_keep_inds,old_grandparents_of_non_root_keep_inds);
            new_weights=1./(1./old_weights_to_parents+1./old_weights_parents_to_grandparents);
            new_weights=[new_weights;new_weights];
        otherwise
            error('Unknown graph reduction method');
    end
    new_W=sparse(i_inds,j_inds,new_weights,new_N,new_N);
    
    % Update parents
    new_root=find(keep_inds==root);
    parents=zeros(size(keep_inds));
    parents([1:new_root-1,new_root+1:end])=new_non_root_parents;
    
    % Update depths
    depths=depths(keep_inds);
    depths=depths/2;
    
    % Store new tree
    Gs{lev+1}=gsp_graph(new_W,Gs{lev}.coords(keep_inds,:),G.limits);
    Gs{lev+1} = gsp_create_laplacian( Gs{lev+1}, G.lap_type);
    Gs{lev+1}.type='tree';
    Gs{lev+1}.root=new_root;
    Gs{lev+1} = gsp_copy_graph_attributes(Gs{lev},0,Gs{lev+1});
    
    if param.compute_full_eigen
        Gs{lev+1}=gsp_compute_fourier_basis(Gs{lev+1});
    end
    
    % Replace current adjacency matrix and root
    old_W=new_W;
    root=new_root;
end


end

  


