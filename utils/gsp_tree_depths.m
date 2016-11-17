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


%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_tree_depths.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.7.0
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

