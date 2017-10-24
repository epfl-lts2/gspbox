function [depths,parents] = gsp_tree_depths(A,root)

if gsp_check_connectivity(A) == 0
    error('Graph is not connected');
end

N = size(A,1);
assigned = root;
depths = zeros(N,1);
parents = zeros(N,1);

next_to_expand = root;
current_depth = 1;

while ( length(assigned) < N )
    new_entries_whole_round = [];
    for i = 1:length(next_to_expand)
        neighbors = find(A(next_to_expand(i),:)>1e-7);
        new_entries = setdiff(neighbors,assigned);
        parents(new_entries) = next_to_expand(i);
        depths(new_entries)=current_depth;
        assigned=[assigned;new_entries'];
        new_entries_whole_round=[new_entries_whole_round;new_entries'];
    end
    current_depth=current_depth+1;
    next_to_expand=new_entries_whole_round;
end

end

