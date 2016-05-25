function [ F ] = gsp_filterbank_matrix(G,g,param )
%GSP_FILTERBANK_MATRIX Create the matrix of the filterbank frame
%   Usage:  F = gsp_filterbank_matrix(G, g param );
%           F = gsp_filterbank_matrix(G,g);
%
%   Input parameters:
%         G     : Graph
%         g     : Filters
%         param : Structure of optional parameter
%   Output parameters:
%         F     : Frame
%
%   This function create the matrix associated to the filterbank g. The
%   size of the matrix is MN x N, where M is the number of filters.
%
%   param a Matlab structure containing the following fields:
% 
%      param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%       By default, it is 1.
%
%
%   Url: http://lts2research.epfl.ch/gsp/doc/filters/gsp_filterbank_matrix.php

% Copyright (C) 2013-2016 Nathanael Perraudin, Johan Paratte, David I Shuman.
% This file is part of GSPbox version 0.6.0
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



%   AUTHOR : Nathanael Perraudin
%   TESTING: test_filter

% Optional input arguments
if nargin<4, param=struct; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end

if param.verbose && G.N>2000
    warning('Create a big matrix, you can use other methods.');
end

Nf = length(g);

Ft = gsp_filter_analysis(G,g,eye(G.N));

F = zeros(size(Ft'));
for ii = 1:Nf
    F(:,G.N*(ii-1)+(1:G.N)) = Ft(G.N*(ii-1)+(1:G.N),:);
end

end


