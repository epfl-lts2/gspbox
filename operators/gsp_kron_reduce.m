function [Greduced]=gsp_kron_reduce(G,keep_inds)
%GSP_KRON_REDUCE Performs Kron reduction
%   Usage: Greduced=gsp_kron_reduce(G,keep_inds);
%          Lreduced=gsp_kron_reduce(L,keep_inds);
%
%   Input parameters:
%         G            : Graph structure or graph Laplacian matrix.
%         keep_inds    : The set of indices to keep in the reduced graph.
%   Output parameters:
%         Greduced     : The Kron-reduced graph structure or Laplacian.
%
%   'gsp_kron_reduce(G,keep_inds)' performs Kron reduction:
%
%   .. L_reduced = L_{V_1,V_1} - L_{V_1,V_2} * [L_{V_2,V_2}]^-1 * L_{V_2,V_1}
%
%   .. math:: {\cal L}_{reduced}={\cal L}_{{\cal V}_1,{\cal V}_1}-{\cal L}_{{\cal V}_1,{\cal V}_2} \left[{\cal L}_{{\cal V}_2,{\cal V}_2}\right]^{-1} {\cal L}_{{\cal V}_2,{\cal V}_1}
%
%   If a matrix is given, then a matrix is returned
%
%   Example:::
%
%           N = 64;
%           param.distribute = 1;
%           param.Nc = 5;
%           param.regular = 1;
%           G = gsp_sensor(N,param);
%           ind = 1:2:N;
%           Gnew = gsp_kron_reduction( G,ind );
%           figure;
%           subplot(121)
%           gsp_plot_graph(G);
%           title('Original graph');
%           subplot(122)
%           gsp_plot_graph(Gnew);
%           title('Kron reduction');
%
%   See also:  
%
%   Notes: may be able to speed this up with LAMG toolbox
%
%   Demos:  
% 
%   References: dorfler2013kron

%   Author : David I Shuman, Nathanael Perraudin.
%   Date   : 26 November 2016
%   Testing: test_operators

if isstruct(G)
    if ~strcmp(G.lap_type,'combinatorial')
        error('Not implemented');
    end
    if G.directed
        error('GSP_KRON_REDUCTION: this only work for undirected graphs.');
    end
    L = G.L;
else
    L = G;
end

N=size(L,1);
if N~= size(L,2)
    error('Graph Laplacian should be a square matrix');
end
elim_inds=setdiff(1:N,keep_inds);

Lreduced=L(keep_inds,keep_inds)-L(keep_inds,elim_inds)*(L(elim_inds,elim_inds)\L(elim_inds,keep_inds));

% Make the laplacian symetric if it is almost symetric!
if sum(sum(abs(Lreduced-transpose(Lreduced))))<eps*sum(sum(abs(Lreduced)))
    Lreduced = (Lreduced+transpose(Lreduced))/2;
end

if isstruct(G)
    Wnew = diag(diag(Lreduced)) - Lreduced;
    Greduced = gsp_copy_graph_attributes(G,0);
    Greduced.W = Wnew;
    Greduced.type='Kron reduction';
    Greduced = gsp_graph_default_parameters(Greduced);
    Greduced.plotting.limits = G.plotting.limits;
    Greduced.coords = G.coords(keep_inds,:);
else % return a Laplacian matrix if a Laplacian matrix was entered
    Greduced=Lreduced;
end

end