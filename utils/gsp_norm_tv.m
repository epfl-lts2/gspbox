function y = gsp_norm_tv(G,x)
%GSP_NORM_TV TV norm on graph
%   Usage:  y = gsp_norm_tv(G,x);
%
%   Input parameters:
%         G     : Graph structure
%         x     : Signal on graph
%   Output parameters:
%         y     : Norm
%
%   Compute the TV norm of a signal on a graph
%
%   See also: gsp_prox_tv
%
%   Url: http://lts2research.epfl.ch/gsp/doc/utils/gsp_norm_tv.php

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

% Author: Nathanael Perraudin
% Date:   25 March 2014
% Testing: test_gsp_prox

if ~isfield(G,'v_in')
    G = gsp_adj2vec(G);
    warning(['GSP_PROX_TV: To be more efficient you should run: ',...
        'G = gsp_adj2vec(G); before using this proximal operator.']);
end


y = sum(abs(gsp_grad(G,x)));

end

