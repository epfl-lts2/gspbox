% function c_est = interpolate_on_complete_graph(c_obs, ind_obs, L, reg, N, solver)
%
% This function takes as inputs: 
% - c_obs : low-dimensional graph signal to interpolate.
% - ind_obs : sampled nodes at which c_obs is observed
% - L is a polynomial in the Laplacian matrix of this graph 
% (can be either a matrix or a function handle g(x) that 
% takes as input a vector of size N and outputs its multiplication
% with the choosen polynomial of L). 
% - N is the number of nodes. 
% - solver is the MATLAB solver to use to solve the optimisation problem. 
% Either 'gmres', or 'cgs'. 
% 
% and outputs:
% - c_est the interpolated signals to the whole graph.
% 
% BEWARE: this code assumes that there is no repetition of indices in ind_obs
% 
% Copyright (C) 2016 Nicolas Tremblay, Gilles Puy.
% This file is part of the CSCbox (Compressive Spectral Clustering toolbox)
%
% The CSCbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The CSCbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% If you use this toolbox please kindly cite
%     N. Tremblay, G. Puy, R. Gribonval and P. Vandergheynst.
%     Compressive Spectral Clustering.
%     ArXiv e-prints:1602.02018 Feb. 2016.

function c_est = interpolate_on_complete_graph(c_obs, ind_obs, L, reg, N, solver)

% Zero-fill c_obs
b = sparse(ind_obs, ones(numel(c_obs), 1), c_obs, N, 1);

% Matrix to invert
MtM = sparse(ind_obs, ind_obs, ones(numel(ind_obs), 1), N, N);
if isa(L, 'function_handle')
    A = @(x) MtM*x+ reg*(L(x));
    % Invert the system
    if strcmp(solver, 'cgs')
        [c_est] = cgs(A, b, 1e-6, 100);
    elseif strcmp(solver, 'gmres')
        [c_est] = gmres(A, b, [],1e-6, 100);
    else
        error(['interpolate_on_complete_graph: solver must be either set to ''cgs'' or to ''gmres''']);
    end
else
    MtM = MtM + reg*L;
    % Invert the system
    if strcmp(solver, 'cgs')
        [c_est] = cgs(MtM, b, 1e-6, 100);
    elseif strcmp(solver, 'gmres')
        [c_est] = gmres(MtM, b, [],1e-6, 100);
    else
        error(['interpolate_on_complete_graph: solver must be either set to ''cgs'' or to ''gmres''']);
    end
end

end
