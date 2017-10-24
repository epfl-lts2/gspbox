function [ sparserL , sparserA, sparserW ] = gsp_graph_sparsify_old(L,epsilon)
% The graph sparsification algorithm of D. A. Spielman and N. Srivastava, 
% "Graph sparsification by effective resistances," SIAM J. Comput., 
% vol. 40, no. 6, pp. 1913-1926, 2011. 

N=size(L,1);

% Epsilon should be between 1/sqrt(N) and 1, with a larger epsilon leading to a sparser graph
if ( (epsilon <= 1/sqrt(N)) || (epsilon >1) )
    error('Epsilon out of required range');
end

% Compute resistance distances, and check if the original graph is
% connected
resistance_distances = gsp_compute_resistance_distances_old(L);
W=diag(diag(L))-L;
W(W<1e-10)=0;
W=sparse(W);
original_connected=gsp_check_connectivity_undirected(W);
if (original_connected==0)
    warning('Original graph not connected before sparsification');
end

% Initialize the probability distribution that will be used to select edges in the sparsified graph
[start_nodes,end_nodes,weights] = find(tril(W));
weights=max(0,weights);
Re=max(0,resistance_distances(sub2ind(size(resistance_distances),start_nodes,end_nodes)));
Pe=weights.*Re;
Pe=Pe/sum(Pe);
    
max_tries=10; % maximum number of tries to get a connected graph

for ii=1:max_tries
    
    % Set Q
    C0=1/30; % Rudelson, 1996 Random Vectors in the Isotropic Position (too hard to figure out actual C0)
    C= 4*C0; % Rudelson and Vershynin, 2007, Thm. 3.1
    q=round(9*C^2*N*log(N)/(epsilon^2));

    % Choose random edges in the new graph according the probability
    % distribution initialized above
    results=gendist(Pe',q,1);
    spin_counts=hist(results,1:length(Pe'));
    per_spin_weights=weights./(q*Pe);

    % Tally the new weights and form the new graph
    new_weights=spin_counts'.*per_spin_weights;
    sparserW=sparse(start_nodes,end_nodes,new_weights,N,N);
    sparserW=sparserW+sparserW';
    sparserL=diag(sum(sparserW))-sparserW;
    sparserD=diag(diag(sparserL));
    sparserW=sparserD-sparserL;
    sparserA=sign(sparserW);

    % Check if new graph is connected. If not, reduce epsilon and try again
    new_graph_connected=gsp_check_connectivity_undirected(sparserA);
    if new_graph_connected
        break;
    elseif ii==max_tries
        warning('Despite attempts to reduce epsilon, sparsified graph is disconnected');
    else
        epsilon=epsilon-(epsilon-1/sqrt(N))/2;
    end
end

end

