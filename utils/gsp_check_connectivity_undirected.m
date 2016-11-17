function [ bool, in ] = gsp_check_connectivity_undirected( G,param )
%GSP_CHECK_CONNECTIVITY_UNDIRECTED Check if the graph G is aperiodic strongly connected
%   Usage: bool = gsp_check_connectivity_undirected( G );
%          bool = gsp_check_connectivity_undirected( L );
%          bool = gsp_check_connectivity_undirected( W );
%          [bool,in]=gsp_check_connectivity_undirected( ... );
%
%   Input parameters:
%       G,W,L: Graph, Laplacian matrix or Weight martrix
%       param: Optional parameters
%
%   Output parameters
%       bool: Boolean
%       in  : Nodes without any in connections
%
%   Test if each node have at least one in connection and one out
%   connection. If this simple test give good results, the function compute
%   the perron vector of G and test it. It might take some time.
%
%   param is an optional structure that contains the following field
%
%    param.verbose*: display parameter - 0 no log - 1 display the errors
%   
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_check_connectivity_undirected.php

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

%TODO: Make this function better

% Date: 6 june 2013
% Author: Nathanael Perraudin

% Handle Input parameters
if nargin<1, error('Not enought inputs parameters!'); end




if nargin<2, param=struct; end
    
if ~isfield(param, 'verbose'), param.verbose = 0; end


if isstruct(G)
    A=G.W;
else
    A=G;
end


% Remove the diagonal
A=A-diag(diag(A));


% Check the connecivity -- simple
bool=~boolean((sum(1.-(sum(A,1)>0))+sum(1.-(sum(A,2)>0))));
in=find(1.-(sum(A,1)>0));
try   
    if bool
        d=bfs(A,1);
        bool=~any(d==-1);
    end
catch
    warning('Using old slow method, install gaimc for speed gain')
    if bool
        D=diag(sum(A).^(-0.5));
        A=D*A*D; % eigenvalues of A between -1 and 1
        % We need to check the number of 1 eigenvalues
        max_eig_P = eigs(speye(size(A))+A,2);
        if abs(max_eig_P(2)-2)<10e-12
            if param.verbose
                fprintf('   ---   the second eigenvalue is 0   ---\n');
            end
            bool=0;
        end
    end

    if param.verbose
        fprintf('   ---   End of tests   ---\n');   
        if bool
            fprintf('    The graph is connected.\n')
        else
            fprintf('    The graph is not connected.\n');
        end
    end
end

end


