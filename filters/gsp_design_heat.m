function g = gsp_design_heat(G,tau,param)
%GSP_DESIGN_HEAT Design a simple heat kernel
%   Usage: gsp_design_heat(G);
%          gsp_design_heat(G,tau);
%          gsp_design_heat(G,tau,param);
%
%   Input parameters:
%       G       : Graph structure
%       tau     : scaling parameter (default 10)
%
%   Output parameters
%       g       : filter
%
%   This function design the following filter:
%
%       g(x) =  exp(-tau*x/lmax) 
%
%   If tau is a vector, the function returns a cell array of filters.
%
%   param is an optional structure containing the following fields
%
%    param.verbose*: verbosity level. 0 no log - 1 display warnings.
%     (default 1) 
%    param.normalize*: Normalize the kernel (works only if the
%     eigenvalues are present in the graph. Use gsp_compute_fourier_basis
%     for this.) (default 0) 
%
%   This function will compute the maximum eigenvalue of the laplacian. To
%   be more efficient, you can precompute it using:
%
%       G = gsp_estimate_lmax(G);
%
%   Example:
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         g = gsp_design_heat(G);   
%         gsp_plot_filter(G,g);  
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_design_heat.php

% Copyright (C) 2013-2014 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.3.0
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


% Author: Nathanael Perraudin
% Date  : 28 July 2014

if nargin < 3
    param = struct;
end

if nargin < 2
    tau = 10;
end

if ~isfield(param,'verbose'), param.verbose = 1; end
if ~isfield(param,'normalize'), param.normalize = 0; end

if numel(tau)>1
    Nf = numel(tau);
    g = cell(Nf,1);
    for ii = 1:Nf
        g{ii} = gsp_design_heat(G,tau(ii),param); 
    end
    return
end

if isstruct(G)
    if ~isfield(G,'lmax')
        if param.verbose
            fprintf('GSP_DESIGN_HEAT has to compute lmax \n')
        end
        G = gsp_estimate_lmax(G);
    end
   lmax = G.lmax;
else
   lmax = G;
end

if param.normalize
    if ~isfield(G,'E')
        error(['GSP_DESIGN_HEAT: You need the eigenvalues',...
            ' to normalize the kernel']);
    end
    gu = @(x) exp(- tau * x/lmax); 
    ng = norm(gu(G.E));
    g = @(x) exp(- tau * x/lmax) / ng; 
else
    % g{1} = @(x) exp(- tau * x/lmax); 
    g = @(x) exp(- tau * x/lmax); 
end

end


