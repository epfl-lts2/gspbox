function [L_reduced]=gsp_kron_reduce_old(L,keep_inds)
%GSP_KRON_REDUCE Kron reduction
%   Usage:  L_reduced=gsp_kron_reduce(L,keep_inds);
%
%   Input parameters:
%         L             : Graph Laplacian.
%         keep_inds     : The set of indices to keep in the reduced graph.
%   Output parameters:
%         L_reduced     : The Kron-reduced graph Laplacian.
%
%   'gsp_kron_reduce(L,keep_inds)' performs Kron reduction:
%
%   .. L_reduced = L_{V_1,V_1} - L_{V_1,V_2} * [L_{V_2,V_2}]^-1 * L_{V_2,V_1}
%
%   .. math:: {\cal L}_{reduced}={\cal L}_{{\cal V}_1,{\cal V}_1}-{\cal L}_{{\cal V}_1,{\cal V}_2} \left[{\cal L}_{{\cal V}_2,{\cal V}_2}\right]^{-1} {\cal L}_{{\cal V}_2,{\cal V}_1}
%
%   where the matrix function $\hat{h}({\cal L})$ is defined as
%
%   See also:  
%
%   Notes: may be able to speed this up with LAMG toolbox
%

%   Demos:  
% 
%   References: F. Doerfler and F. Bullo, "Kron reduction of graphs with
%   applications to electrical networks," IEEE Trans. Circuits and
%   Systems-I: Regular Papers, vol. 60, no. 1 pp. 150-163, Jan. 2013.

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE: 

N=size(L,1);
if N~= size(L,2)
    error('Graph Laplacian should be a square matrix');
end
elim_inds=setdiff(1:N,keep_inds);
L_reduced=L(keep_inds,keep_inds)-L(keep_inds,elim_inds)*(L(elim_inds,elim_inds)\L(elim_inds,keep_inds));